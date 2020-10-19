// Reaktoro is a unified framework for modeling chemically reactive systems.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    double t0 = 0;
    double dt = 0.1; // in seconds
    int n = 1;
    Index ncells = 1;

    double initial_sa_m2mol = 1.0; // defined in PHREEQC, m2/mol
    double initial_sa_m2kg;

    double molar_mass_halite = 58.44; // g / mol
    double molar_mass_calcite = 100.0869; // g / mol
    double molar_mass_dolomite = 184.40; // g / mol
    double molar_mass_quartz = 60.083; // g / mol
    double molar_mass_kfeldspar = 278.3315; // g / mol
    double molar_mass_kaolinite = 258.071; // g / mol // or 258.1604

    double initial_mineral_amount_mol = 10.0;

    //std::string activity_model = "pitzer";
    std::string activity_model = "hkf";
    std::cout << "Activity model: " << activity_model << std::endl;

    std::string result_filename = "kinetics-simple-benchmark-6-kin-0-eq-dt-" + std::to_string(dt) +
                                                                        "-n-" +  std::to_string(n) +
                                                                        "-ncells-" +  std::to_string(ncells) +
                                                                        "-activity-model-" + activity_model +
                                                                        "-uniform-sa-1.0-m2mol.txt";

    Database database("supcrt07.xml");

    DebyeHuckelParams dhModel{};
    dhModel.setPHREEQC();

    ChemicalEditor editor(database);
    if(activity_model=="pitzer")
        editor.addAqueousPhaseWithElements("H Cl S O Ca Na K Mg C Si Al")
                .setChemicalModelDebyeHuckel(dhModel)
                .setChemicalModelPitzerHMW()
                .setActivityModelDrummondCO2();
    else if(activity_model=="hkf")
        editor.addAqueousPhaseWithElements("H Cl S O Ca Na K Mg C Si Al")
                .setChemicalModelDebyeHuckel(dhModel);

    std::vector<std::string> minerals = {"Halite", "Calcite", "Dolomite", "K-Feldspar", "Quartz", "Kaolinite"};
    for(const auto& mineral : minerals)
        editor.addMineralPhase(mineral);

    ChemicalSystem system(editor);

    const auto I = ChemicalProperty::ionicStrength(system);

    /*
     * # Source: Palandri & Kharaka (2004)
        # The following parameters are specified in KINETICS for each kinetic mineral in the system
        #    PARM(1): initial reactive surface area (m2/mole) - History matching parameter if M0 > 0
        #    PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        #    PARM(3): used when M0 = 0, Nucleation surface area
        #    PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        #    PARM(5): when M0 > 0, rate constant ratio (K_precipitation / kappa_dissolution), (=1 if no data)

        # k25(a,n,b): Rate constant k at 25oC, pH = 0, mole m-2 s-1.
        # Ea(a,n,b):  Arrhenius activation energy E, J mole-1.
        # n(a,b):     Reaction order n with respect to H+ (apart from carbonates and sulfides).
        #           (a) = acid mechanism (H+)
        #           (n) = neutral mechanism
        #           (b) = base mechanism (OH-), apart from carbonates (HCO3-)
        # Note: for sulfide, other mechanisms (Fe+3 and O2)
        #       for carbonates, followed Xu et al (2010)
     */
    initial_sa_m2kg = initial_sa_m2mol / molar_mass_halite / 1e-3;
    std::string eq_str_halite = "Halite = Na+ + Cl-";
    MineralReaction min_reaction_halite = editor.addMineralReaction("Halite")
            .setEquation(eq_str_halite)
            .addMechanism("logk = -0.21 mol/(m2*s); Ea = 7.4 kJ/mol")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg"); // conversion from
    Reaction reaction_halite = createReaction(min_reaction_halite, system);
    reaction_halite.setName("Halite reaction");
    ReactionRateFunction rate_func_halite_llnl = [&min_reaction_halite, &reaction_halite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_halite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_halite.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_halite.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        double Omega = exp(lnOmega).val;
        double Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();

        // The initial composition of the chemical system
        const auto n0_rxn = reaction_halite.initialAmounts();
//        std::cout << "n0_rxn : " << tr(n0_rxn) << std::endl;
//        getchar();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Cl")
        const Index i_cl = system.indexElementWithError("Cl");
        const auto tot_cl = b(i_cl);

        // Calculate TOT("Na")
        const Index i_na = system.indexElementWithError("Na");
        const auto tot_na = b(i_na);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_halite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_halite.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // k25a =   0
        // Eaa =    0
        // na = 0
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.0;
        // k25n =   0.616595002
        // Ean =    7400
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 0.616595002 * exp(- 7.4 / R * (1.0/T - 1.0/298.15));
        // k25b =   0
        // Eab =    0
        // nb = 0
        const auto kappa_base = 0.0;

        // Precipitation and dissolution constant rates
        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;

        /*
         * Halite
        -start
        10 k25a =	0
        20 k25n =	0.616595002
        30 k25b =	0
        40 Eaa =	0
        50 Ean =	7400
        60 Eab =	0
        70 na =	0
        80 nb =	0
        90  SatIndex = SI("Halite")
        100 SatRatio = SR("Halite")
        110 R = 8.314472 # (gas constant J/mol/Kelvin)
        120 REM If SI too small, then No dissolution nor precipitation
        130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        135 IF (SatRatio > 10000) THEN SatRatio = 10000
        140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        160 t = 1    # surface area correction factor
        170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        200 REM acid mechanism
        210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        220 REM neutral mechanism
        230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        240 REM base mechanism
        250   IF (k25b > 0) THEN KTb = k25b * EXP((-Eab/R) * (1/TK - 1/298.15)) * ACT("H+")^nb
        300 REM Dissolution
        310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
        320 IF (SatIndex > 0) THEN GOTO 400
        330   REM Do not dissolve more than what is available
        340   IF (moles > M) THEN moles = M
        350   GOTO 900
        400 REM Precipitation when M0 > 0
        410 IF (M0 = 0) THEN GOTO 500
        420   moles = moles * PARM(5)
        430   GOTO 600
        500 REM Precipitation when M0 = 0 (secondary mineral)
        520   surfArea = 1.0e-5  # default nucleation surface area
        530   IF (M > 0) THEN surfArea = M * PARM(3)
        540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
        600 REM Do not precipitate more than the elements in solution
        610   maxMol = TOT("Na")
        620   IF (maxMol > TOT("Cl")) THEN maxMol =TOT("Cl")
        690   IF (maxMol < -moles) THEN moles = -maxMol
        900 SAVE moles
        -end
         */

        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res +=  nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

            // REM Do not precipitate more than the elements in solution
            // maxMol = TOT("Na")
            // IF (maxMol > TOT("Cl")) THEN maxMol =TOT("Cl")
            // IF (maxMol < -moles) THEN moles = -maxMol
            // SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_na;
//            if(max_mol > tot_cl) max_mol = tot_cl;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);

        }

        return res;

    };
    reaction_halite.setRate(rate_func_halite_llnl);

    initial_sa_m2kg = initial_sa_m2mol / molar_mass_calcite / 1e-3;
    std::string eq_str_calcite = "Calcite + H+  =  Ca++ + HCO3-";
    MineralReaction min_reaction_calcite = editor.addMineralReaction("Calcite")
            .setEquation(eq_str_calcite)
            .addMechanism("logk = -5.81 mol/(m2*s); Ea = 23.5 kJ/mol")
            .addMechanism("logk = -0.30 mol/(m2*s); Ea = 14.4 kJ/mol; a[H+] = 1.0")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");
    Reaction reaction_calcite = createReaction(min_reaction_calcite, system);
    reaction_calcite.setName("Calcite reaction");
    ReactionRateFunction rate_func_calcite_llnl = [&min_reaction_calcite, &reaction_calcite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_calcite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_calcite.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_calcite.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        double Omega = exp(lnOmega).val;
        double Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_calcite.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Ca")
        const Index i_ca = system.indexElementWithError("Ca");
        const auto tot_ca = b(i_ca);

        // Calculate TOT("C")
        const Index i_c = system.indexElementWithError("C");
        const auto tot_c = b(i_c);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_calcite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        double ionic_strength = I(properties).val;

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_calcite.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // k25a =	0.501187234
        // Eaa =	14400
        // na =	1
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.501187234 * exp(- 14.4 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 1.0);
        // k25n =	1.54882E-06
        // Ean =	23500
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 1.54882e-06 * exp(- 23.5 / R * (1.0/T - 1.0/298.15));
        // k25b =	0
        // Eab =	0
        // nb =	0
        const auto kappa_base = 0.0;

        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;

        /*
         * Calcite
            -start
            10 k25a =	0.501187234
            20 k25n =	1.54882E-06
            30 k25b =	0
            40 Eaa =	14400
            50 Ean =	23500
            60 Eab =	0
            70 na =	1
            80 nb =	0
            90  SatIndex = SI("Calcite")
            100 SatRatio = SR("Calcite")
            110 R = 8.314472 # (gas constant J/mol/Kelvin)
            120 REM If SI too small, then No dissolution nor precipitation
            130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
            135 IF (SatRatio > 10000) THEN SatRatio = 10000
            140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
            160 t = 1    # surface area correction factor
            170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
            200 REM acid mechanism
            210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
            220 REM neutral mechanism
            230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
            240 REM base mechanism
            250   KTb = 0     # Xu(2010)
            300 REM Dissolution
            310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            320 IF (SatIndex > 0) THEN GOTO 400
            330   REM Do not dissolve more than what is available
            340   IF (moles > M) THEN moles = M
            350   GOTO 900
            400 REM Precipitation when M0 > 0
            410 IF (M0 = 0) THEN GOTO 500
            420   moles = moles * PARM(5)
            430   GOTO 600
            500 REM Precipitation when M0 = 0 (secondary mineral)
            520   surfArea = 1.0e-5  # default nucleation surface area
            530   IF (M > 0) THEN surfArea = M * PARM(3)
            540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
            600 REM Do not precipitate more than the elements in solution
            610   maxMol = TOT("Ca")
            620   IF (maxMol > TOT("C")) THEN maxMol =TOT("C")
            690   IF (maxMol < -moles) THEN moles = -maxMol
            900 SAVE moles
            -end
        -end
         */

        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res +=  nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

//            600 REM Do not precipitate more than the elements in solution
//            610   maxMol = TOT("Ca")
//            620   IF (maxMol > TOT("C")) THEN maxMol =TOT("C")
//            690   IF (maxMol < -moles) THEN moles = -maxMol
//            900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_ca;
//            if(max_mol > tot_c) max_mol = tot_c;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);

        }

    };
    reaction_calcite.setRate(rate_func_calcite_llnl);

    initial_sa_m2kg = initial_sa_m2mol / molar_mass_dolomite / 1e-3;
    std::string eq_str_dolomite = "Dolomite  + 2*H+  =  Ca++ + Mg++ + 2*HCO3-";
    MineralReaction min_reaction_dolomite = editor.addMineralReaction("Dolomite")
            .setEquation(eq_str_dolomite)
            .addMechanism("logk = -7.53 mol/(m2*s); Ea = 52.2 kJ/mol")
            .addMechanism("logk = -3.19 mol/(m2*s); Ea = 36.1 kJ/mol; a[H+] = 0.5")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");
    Reaction reaction_dolomite = createReaction(min_reaction_dolomite, system);
    reaction_dolomite.setName("Dolomite reaction");
    ReactionRateFunction rate_func_dolomite_llnl = [&min_reaction_dolomite, &reaction_dolomite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_dolomite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_dolomite.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_dolomite.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        auto Omega = exp(lnOmega).val;
        auto Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_dolomite.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Mg")
        const Index i_mg = system.indexElementWithError("Mg");
        const auto tot_mg = b(i_mg);

        // Calculate TOT("C")
        const Index i_c = system.indexElementWithError("C");
        const auto tot_c = b(i_c);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_dolomite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_dolomite.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // k25n =	2.95121E-08
        // Ean =	52200
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 2.95121E-08 * exp(- 52.2 / R * (1.0/T - 1.0/298.15));
        // k25a =	0.000645654
        // Eaa =	36100
        // na =	0.5
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.000645654 * exp(- 36.1 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.5);
        // k25b =	0
        // Eab =	0
        // nb =	0
        const auto kappa_base = 0.0;

        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;
        /*
         * Dolomite-dis
        -start
        10 k25a =	0.000645654
        20 k25n =	2.95121E-08
        30 k25b =	0
        40 Eaa =	36100
        50 Ean =	52200
        60 Eab =	0
        70 na =	0.5
        80 nb =	0
        90  SatIndex = SI("Dolomite-dis")
        100 SatRatio = SR("Dolomite-dis")
        110 R = 8.314472 # (gas constant J/mol/Kelvin)
        120 REM If SI too small, then No dissolution nor precipitation
        130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        135 IF (SatRatio > 10000) THEN SatRatio = 10000
        140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        160 t = 1    # surface area correction factor
        170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        200 REM acid mechanism
        210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        220 REM neutral mechanism
        230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        240 REM base mechanism
        250   KTb = 0     # Xu(2010)
        300 REM Dissolution
        310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
        320 IF (SatIndex > 0) THEN GOTO 400
        330   REM Do not dissolve more than what is available
        340   IF (moles > M) THEN moles = M
        350   GOTO 900
        400 REM Precipitation when M0 > 0
        410 IF (M0 = 0) THEN GOTO 500
        420   moles = moles * PARM(5)
        430   GOTO 600
        500 REM Precipitation when M0 = 0 (secondary mineral)
        520   surfArea = 1.0e-5  # default nucleation surface area
        530   IF (M > 0) THEN surfArea = M * PARM(3)
        540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
        600 REM Do not precipitate more than the elements in solution
        610   maxMol = TOT("Ca")
        620   IF (maxMol > TOT("C")/2) THEN maxMol =TOT("C")/2
        630   IF (maxMol > TOT("Mg")) THEN maxMol =TOT("Mg")
        690   IF (maxMol < -moles) THEN moles = -maxMol
        900 SAVE moles
        -end
         */

        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res += nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

//            600 REM Do not precipitate more than the elements in solution
//            610   maxMol = TOT("Ca")
//            620   IF (maxMol > TOT("C")/2) THEN maxMol =TOT("C")/2
//            630   IF (maxMol > TOT("Mg")) THEN maxMol =TOT("Mg")
//            690   IF (maxMol < -moles) THEN moles = -maxMol

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_ca;
//            if(max_mol > tot_c/2) max_mol = tot_c /2;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);

        }

        return res;

    };
    reaction_dolomite.setRate(rate_func_dolomite_llnl);

    // K-Feldspar: K(AlSi3)O8 # potassium feldspar, orthoclase
    initial_sa_m2kg = initial_sa_m2mol / molar_mass_kfeldspar / 1e-3;
    // KAlSi3O8 + 4 H+  =  Al+++ + K+ + 2 H2O + 3 SiO2
    std::string eq_str_kfeldspar = "K-Feldspar + 4*H+ = Al+++ + 3*SiO2(aq) + K+ + 2*H2O(l)";
    MineralReaction min_reaction_kfeldspar = editor.addMineralReaction("K-Feldspar")
            .setEquation("K-Feldspar + H2O(l) + H+ = Al+++ + 3*HSiO3- + K+")
            .addMechanism("logk = -12.41 mol/(m2*s); Ea = 38.0 kJ/mol")
            .addMechanism("logk = -10.06 mol/(m2*s); Ea = 51.7 kJ/mol; a[H+] = 0.5")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");
    Reaction reaction_kfeldspar = createReaction(min_reaction_kfeldspar, system);
    reaction_kfeldspar.setName("K-Feldspar reaction");
    ReactionRateFunction rate_func_kfeldspar_llnl = [&min_reaction_kfeldspar, &reaction_kfeldspar, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_kfeldspar.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_kfeldspar.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_kfeldspar.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        auto Omega = exp(lnOmega).val;
        auto Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_kfeldspar.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        const Index i_si = system.indexElementWithError("Si");
        const auto tot_si = b(i_si);

        // Calculate TOT("K")
        const Index i_k = system.indexElementWithError("K");
        const auto tot_k = b(i_k);

        // Calculate TOT("Al")
        const Index i_al = system.indexElementWithError("Al");
        const auto tot_al = b(i_al);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_kfeldspar.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_kfeldspar.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // Acid mechanism
        // k25a =  2.51189E-11
        // Eaa =    57100
        // na = 0.5
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 2.51189E-11 * exp(- 57.1 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, 0.5);

        // Neutral mechanism
        // k25n =   3.89045E-13
        // Ean =    38000
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 3.89045E-13 * exp(- 38.0 / R * (1.0/T - 1.0/298.15));

        // Base mechanism
        // k25b =   6.30957E-22
        // Eab =    94100
        // nb = -0.823
        const auto kappa_base = 6.30957E-22 * exp(- 94.1 / R * (1.0/T - 1.0/298.15)) * std::pow(activity_h, -0.823);

        // Define dissolution and precipitation rates
        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;

        /*
         * K-Feldspar
            -start
            10 k25a =	2.51189E-11
            20 k25n =	3.89045E-13
            30 k25b =	6.30957E-22
            40 Eaa =	57100
            50 Ean =	38000
            60 Eab =	94100
            70 na =	0.5
            80 nb =	-0.823
            90  SatIndex = SI("K-Feldspar")
            100 SatRatio = SR("K-Feldspar")
            110 R = 8.314472 # (gas constant J/mol/Kelvin)
            120 REM If SI too small, then No dissolution nor precipitation
            130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
            135 IF (SatRatio > 10000) THEN SatRatio = 10000
            140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
            160 t = 1    # surface area correction factor
            170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
            200 REM acid mechanism
            210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
            220 REM neutral mechanism
            230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
            240 REM base mechanism
            250   IF (k25b > 0) THEN KTb = k25b * EXP((-Eab/R) * (1/TK - 1/298.15)) * ACT("H+")^nb
            300 REM Dissolution
            310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            320 IF (SatIndex > 0) THEN GOTO 400
            330   REM Do not dissolve more than what is available
            340   IF (moles > M) THEN moles = M
            350   GOTO 900
            400 REM Precipitation when M0 > 0
            410 IF (M0 = 0) THEN GOTO 500
            420   moles = moles * PARM(5)
            430   GOTO 600
            500 REM Precipitation when M0 = 0 (secondary mineral)
            520   surfArea = 1.0e-5  # default nucleation surface area
            530   IF (M > 0) THEN surfArea = M * PARM(3)
            540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
            600 REM Do not precipitate more than the elements in solution
            610   maxMol = TOT("Si")/3
            620   IF (maxMol > TOT("Al")) THEN maxMol =TOT("Al")
            640   IF (maxMol > TOT("K")) THEN maxMol =TOT("K")
            690   IF (maxMol < -moles) THEN moles = -maxMol
            900 SAVE moles
            -end
         */

        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res += nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

            // 600 REM Do not precipitate more than the elements in solution
            // 610   maxMol = TOT("Si")/3
            // 620   IF (maxMol > TOT("Al")) THEN maxMol =TOT("Al")
            // 640   IF (maxMol > TOT("K")) THEN maxMol =TOT("K")
            // 690   IF (maxMol < -moles) THEN moles = -maxMol
            // 900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_si;
//            if(max_mol > tot_al) max_mol = tot_al;
//            if(max_mol > tot_k) max_mol = tot_k;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0)

        }

        return res;

    };
    reaction_kfeldspar.setRate(rate_func_kfeldspar_llnl);

    // Kaolinite: Al2Si2O5(OH)4
    initial_sa_m2kg = initial_sa_m2mol / molar_mass_kaolinite / 1e-3;
    //std::string eq_str_kaolinite = "Kaolinite + 4*H+ = 3*H2O(l) + 2*HSiO3- + 2*Al+++";
    // Al2Si2O5(OH)4 + 6 H+  =  2 Al+++ + 2 SiO2 + 5 H2O
    std::string eq_str_kaolinite = "Kaolinite + 6*H+ = 2*Al+++ + 2*SiO2(aq) + 5*H2O(l)";
    MineralReaction min_reaction_kaolinite = editor.addMineralReaction("Kaolinite")
            .setEquation(eq_str_kaolinite)
            .addMechanism("logk = -13.18 mol/(m2*s); Ea = 22.2 kJ/mol")
            .addMechanism("logk = -11.31 mol/(m2*s); Ea = 65.9 kJ/mol; a[H+] = 0.777")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");
    Reaction reaction_kaolinite = createReaction(min_reaction_kaolinite, system);
    reaction_kaolinite.setName("Kaolinite reaction");
    ReactionRateFunction rate_func_kaolinite_llnl = [&min_reaction_kaolinite, &reaction_kaolinite, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_kaolinite.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_kaolinite.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_kaolinite.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        auto Omega = exp(lnOmega).val;
        auto Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_kaolinite.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        const Index i_si = system.indexElementWithError("Si");
        const auto tot_si = b(i_si);

        // Calculate TOT("Al")
        const Index i_al = system.indexElementWithError("Al");
        const auto tot_al = b(i_al);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_kaolinite.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");
        double activity_h = std::exp(lna(i_h));
        double ionic_strength = I(properties).val;

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_kaolinite.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // Acid mechanism
        // k25a =  0
        // Eaa =   0
        // na = 0
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.0;

        // Neutral mechanism
        // k25n =   1E-13
        // Ean =    7100
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 1E-13 * exp(- 7.1 / R * (1.0/T - 1.0/298.15));

        // Base mechanism
        // k25b =   0
        // Eab =    0
        // nb = 0
        const auto kappa_base = 0.0;

        // Define dissolution and precipitation rates
        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;
        /*
         * Kaolinite
            -start
            10 k25a =	0
            20 k25n =	1E-13
            30 k25b =	0
            40 Eaa =	0
            50 Ean =	7100
            60 Eab =	0
            70 na =	0
            80 nb =	0
            90  SatIndex = SI("Kaolinite")
            100 SatRatio = SR("Kaolinite")
            110 R = 8.314472 # (gas constant J/mol/Kelvin)
            120 REM If SI too small, then No dissolution nor precipitation
            130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
            135 IF (SatRatio > 10000) THEN SatRatio = 10000
            140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
            160 t = 1    # surface area correction factor
            170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
            200 REM acid mechanism
            210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
            220 REM neutral mechanism
            230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
            240 REM base mechanism
            250   IF (k25b > 0) THEN KTb = k25b * EXP((-Eab/R) * (1/TK - 1/298.15)) * ACT("H+")^nb
            300 REM Dissolution
            310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            320 IF (SatIndex > 0) THEN GOTO 400
            330   REM Do not dissolve more than what is available
            340   IF (moles > M) THEN moles = M
            350   GOTO 900
            400 REM Precipitation when M0 > 0
            410 IF (M0 = 0) THEN GOTO 500
            420   moles = moles * PARM(5)
            430   GOTO 600
            500 REM Precipitation when M0 = 0 (secondary mineral)
            520   surfArea = 1.0e-5  # default nucleation surface area
            530   IF (M > 0) THEN surfArea = M * PARM(3)
            540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
            600 REM Do not precipitate more than the elements in solution
            610   maxMol = TOT("Si")/2
            620   IF (maxMol > TOT("Al")/2) THEN maxMol =TOT("Al")/2
            690   IF (maxMol < -moles) THEN moles = -maxMol
            900 SAVE moles
            -end
         */
        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res += nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

            // 600 REM Do not precipitate more than the elements in solution
            // 610   maxMol = TOT("Si")/2
            // 620   IF (maxMol > TOT("Al")/2) THEN maxMol =TOT("Al")/2
            // 690   IF (maxMol < -moles) THEN moles = -maxMol
            // 900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_si;
//            if(max_mol > tot_al) max_mol = tot_al;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0)

        }

        return res;

    };
    reaction_kaolinite.setRate(rate_func_kaolinite_llnl);

    initial_sa_m2kg = initial_sa_m2mol / molar_mass_quartz / 1e-3;
    //std::string eq_str_quartz = "Quartz + H2O(l) = H+ + HSiO3-";
    std::string eq_str_quartz = "Quartz = SiO2(aq)";
    MineralReaction min_reaction_quartz = editor.addMineralReaction("Quartz")
            .setEquation(eq_str_quartz)
            .addMechanism("logk = -13.99 mol/(m2*s); Ea = 87.7 kJ/mol")
            .setSpecificSurfaceArea(initial_sa_m2kg, "m2/kg");
    Reaction reaction_quartz = createReaction(min_reaction_quartz, system);
    reaction_quartz.setName("Quartz reaction");
    ReactionRateFunction rate_func_quartz_llnl = [&min_reaction_quartz, &reaction_quartz, &system, &I](const ChemicalProperties& properties) -> ChemicalScalar {

        // The number of chemical species in the system
        const unsigned num_species = system.numSpecies();

        // The mineral reaction rate
        ChemicalScalar res(num_species, 0.0);
        // Auxiliary variable for calculating mineral reaction rate
        ChemicalScalar f(num_species, 1.0);

        // The universal gas constant (in units of kJ/(mol*K))
        const double R = 8.3144621e-3;

        // The temperature and pressure of the system
        const Temperature T = properties.temperature();

        // Create a Reaction instance
        Reaction reaction(min_reaction_quartz.equation(), system);

        // Calculate the saturation index of the mineral
        const auto lnK = reaction.lnEquilibriumConstant(properties);
        const auto lnQ = reaction.lnReactionQuotient(properties);
        const auto lnOmega = lnQ - lnK;

        const auto lnK_rxn = reaction_quartz.lnEquilibriumConstant(properties);
        const auto lnQ_rxn = reaction_quartz.lnReactionQuotient(properties);
        const auto lnOmega_rxn = lnQ_rxn - lnK_rxn;

        // Calculate the saturation index
        //auto Omega = exp(lnOmega);
        //auto Omega_rxn = exp(lnOmega_rxn);
        auto Omega = exp(lnOmega).val;
        auto Omega_rxn = exp(lnOmega_rxn).val;

        // The composition of the chemical system
        const auto n = properties.composition();
        // The initial composition of the chemical system
        const auto n0_rxn = reaction_quartz.initialAmounts();

        // Amount of elements
        const auto b = system.elementAmounts(n.val);

        // Calculate TOT("Si")
        const Index i_si = system.indexElementWithError("Si");
        const auto tot_si = b(i_si);

        // The index of the mineral
        const Index imineral = system.indexSpeciesWithError(min_reaction_quartz.mineral());

        // The current and the initial number of moles of the mineral
        auto nm = n[imineral];
        auto nm0_rxn = n0_rxn[imineral];

        VectorConstRef lna = properties.lnActivities().val;
        const Index i_h = system.indexSpeciesWithError("H+");

        // The molar mass of the mineral (in units of kg/mol)
        const double molar_mass = system.species(imineral).molarMass();

        // Specific surface area
        const auto ssa = min_reaction_quartz.specificSurfaceArea();

        // Specific surface area correction factor
        auto ssa_corr_factor = 1.0;

        // REM If SI too small, then No dissolution nor precipitation
        // IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        if(lnOmega > -1e-15 && lnOmega < 1e-15)
            return res;
        //  IF (SatIndex > 0) THEN GOTO 900	# IDISPRE=1 in TR model (i.e. only dissolution allowed)
        if(lnOmega > 0)
            return res;
        // IF (SatRatio > 10000) THEN SatRatio = 10000
        if(Omega > 1e4)
            Omega = 1e4;
        // IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        // PARM(2): used when M0 = 0, (Nucleation) SatIndex threshold for precipitation
        double sat_index_threshold = 1e6;
        if(nm.val == 0 && lnOmega < sat_index_threshold)
            return res;
        // IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        // PARM(4): exponent for M/M0 for surface area correction (=1 if no data)
        if(nm0_rxn > 0 && nm.val > 0)
            ssa_corr_factor = pow(nm.val/nm0_rxn, 2/3);

        // k25a =   0
        // Eaa =    0
        // na = 0
        // IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        const auto kappa_acid = 0.0;
        // k25n =	4.52E-14
        // Ean =	90100
        // IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        const auto kappa_neu = 4.52E-14 * exp(- 90.1 / R * (1.0/T - 1.0/298.15));
        // k25b =   0
        // Eab =    0
        // nb = 0
        const auto kappa_base = 0.0;

        // Precipitation and dissolution constant rates
        auto kappa_diss = kappa_neu + kappa_acid + kappa_base;
        auto kappa_pre = kappa_neu + kappa_acid + kappa_base;

        /*
         * Quartz
        -start
        10 k25a =	0
        20 k25n =	4.52E-14 # NB we have not implemented the rate correction below pH=3 and pH=8 (mineral record 3 in TR model), since pH stays within this interval
        30 k25b =	0
        40 Eaa =	0
        50 Ean =	90100
        60 Eab =	0
        70 na =	0
        80 nb =	0
        90  SatIndex = SI("Quartz")
        100 SatRatio = SR("Quartz")
        110 R = 8.314472 # (gas constant J/mol/Kelvin)
        120 REM If SI too small, then No dissolution nor precipitation
        130 IF (SatIndex > -1e-15) AND ( SatIndex < 1e-15) THEN GOTO 900
        132 IF (SatIndex > 0) THEN GOTO 900	# IDISPRE=1 in TR model (i.e. only dissolution allowed)
        135 IF (SatRatio > 10000) THEN SatRatio = 10000
        140 IF (M = 0) AND (SatIndex < PARM(2)) THEN GOTO 900
        160 t = 1    # surface area correction factor
        170 IF (M0 > 0) AND (M > 0) THEN t = (M / M0)^PARM(4)
        200 REM acid mechanism
        210   IF (k25a > 0) THEN KTa = k25a * EXP((-Eaa/R) * (1/TK - 1/298.15)) * ACT("H+")^na
        220 REM neutral mechanism
        230   IF (k25n>0) THEN KTn = k25n * EXP((-Ean/R) * (1/TK - 1/298.15))
        240 REM base mechanism
        250   IF (k25b > 0) THEN KTb = k25b * EXP((-Eab/R) * (1/TK - 1/298.15)) * ACT("H+")^nb
        300 REM Dissolution
        310   moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
        320 IF (SatIndex > 0) THEN GOTO 400
        330   REM Do not dissolve more than what is available
        340   IF (moles > M) THEN moles = M
        350   GOTO 900
        400 REM Precipitation when M0 > 0
        410 IF (M0 = 0) THEN GOTO 500
        420   moles = moles * PARM(5)
        430   GOTO 600
        500 REM Precipitation when M0 = 0 (secondary mineral)
        520   surfArea = 1.0e-5  # default nucleation surface area
        530   IF (M > 0) THEN surfArea = M * PARM(3)
        540   moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
        600 REM Do not precipitate more than the elements in solution
        610   maxMol = TOT("Si")
        690   IF (maxMol < -moles) THEN moles = -maxMol
        900 SAVE moles
        -end
         */

        if(Omega <= 1) // dissolution kinetics
        {
            // moles = PARM(1) * M0 * t * (KTa+KTn+KTb) * (1 - SatRatio) * TIME
            res += f * ssa * nm0_rxn * ssa_corr_factor * molar_mass * kappa_pre * (1 - Omega);

            // Do not dissolve more than what is available
            // 320 IF (SatIndex > 0) THEN GOTO 400
            // 330   REM Do not dissolve more than what is available
            // 340   IF (moles > M) THEN moles = M
            // 350   GOTO 900
//            double total_moles = nm.val; // current amount of mols of available minerals
//            if (lnOmega <= 0 && res > nm.val)
//                res += nm.val;
        }
        else // precipitation kinetics
        {
            // REM Precipitation when M0 = 0 (secondary mineral)
            if(nm0_rxn != 0)
            {
                // Default nucleation surface area
                auto nucl_sa = 1.0e-5;
                // IF (M > 0) THEN surfArea = M * PARM(3)
                if(nm.val) nucl_sa *= nm.val;

                // moles = surfArea * (KTa+KTn+KTb)* PARM(5) * (1 - SatRatio) * TIME
                res += f * nucl_sa * kappa_pre * (1 - Omega);
            }
            else // IF (M0 = 0) THEN GOTO 500
                // Set nucleation rate
                // moles = moles * PARM(5)
                res *= f * kappa_pre / kappa_diss;

            //     600 REM Do not precipitate more than the elements in solution
            // 610   maxMol = TOT("Si")
            // 690   IF (maxMol < -moles) THEN moles = -maxMol
            // 900 SAVE moles

            // TODO: implement if in the kinetic solver
//            // Implement upper bound in precipitation kinetics
//            auto max_mol = tot_si;
//            if(max_mol < -nm.val) return ChemicalScalar(num_species, 0.0);

        }

        return res;

    };
    reaction_quartz.setRate(rate_func_quartz_llnl);

    //ReactionSystem reactions(editor);
    std::vector<Reaction> reaction_minerals = {reaction_halite, reaction_calcite, reaction_dolomite, reaction_quartz, reaction_kfeldspar, reaction_kaolinite};
    ReactionSystem reactions(system, reaction_minerals);

    Partition partition(system);
    partition.setKineticPhases({"Halite", "Calcite", "Dolomite", "K-Feldspar", "Quartz", "Kaolinite"});

    double T = 61.0 + 273.15;       // temperature (in units of celsius)
    double P = 1 * 1.01325 * 1e5;   // pressure (in units of atm)
    EquilibriumInverseProblem problem(partition);
    problem.setTemperature(T);
    problem.setPressure(P);
    problem.add("H2O", 1, "kg");
    problem.add("Na", 6.27 , "mol"); // 6.27 mol / kgw
    problem.add("Cl", 6.27, "mol");  // 6.27 mol / kgw
    problem.add("Mg", 1e-4, "mol");  // 0.0001 mol / kgw
    problem.add("Ca", 1e-4, "mol");  // 0.0001 mol / kgw
    problem.add("C", 1e-4, "mol");  // 0.0001 mol / kgw
    problem.add("Si", 1e-4, "mol");  // 0.0001 mol / kgw
    problem.add("Al", 1e-4, "mol");  // 6.27 mol / kgw
    problem.add("K", 1e-9, "mol");

    problem.pH(7.0);

    //problem.pE(4.0); // the result of PHREEQC has pE = 4.0

    EquilibriumInverseSolver solver(partition);
    ChemicalState state0(system);
    EquilibriumResult result = solver.solve(state0, problem);

    //std::cout << "state0 = \n" << state0 << std::endl;

    ChemicalProperties props = state0.properties();
    auto evaluate_pH = ChemicalProperty::pH(system);
    auto evaluate_alkalinity = ChemicalProperty::alkalinity(system);
    auto pH = evaluate_pH(props);
    auto alkalinity = evaluate_alkalinity(props);

    std::cout << "# iter = " << result.optimum.iterations << std::endl;
    std::cout << "pH                 = " << pH.val << std::endl;
    std::cout << "Alkalinity [eq/L]  = " << alkalinity << std::endl;

    /*
    state0 =
    =================================================================================================================================================================================================================================
    Temperature [K]          Temperature [C]          Pressure [Pa]            Pressure [bar]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    334.15                   61                       101325                   1.01325
    =================================================================================================================================================================================================================================
    Element                  Amount [mol]             Aqueous [mol]            Halite [mol]             Calcite [mol]            Dolomite [mol]           K-Feldspar [mol]         Quartz [mol]             Kaolinite [mol]          Dual Potential [kJ/mol]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Al                       0.0001                   0.0001                   0                        0                        0                        1e-10                    0                        2e-10                    -430.233
    C                        0.0001                   0.0001                   0                        1e-10                    2e-10                    0                        0                        0                        49.173
    Ca                       0.0001                   0.0001                   0                        1e-10                    1e-10                    0                        0                        0                        -498.783
    Cl                       6.27                     6.27                     1e-10                    0                        0                        0                        0                        0                        -172.304
    H                        111.017                  111.017                  0                        0                        0                        0                        0                        4e-10                    -2.35331
    K                        1.1e-09                  1e-09                    0                        0                        0                        1e-10                    0                        0                        -303.023
    Mg                       0.0001                   0.0001                   0                        0                        1e-10                    0                        0                        0                        -396.675
    Na                       6.27                     6.27                     1e-10                    0                        0                        0                        0                        0                        -218.465
    O                        55.5084                  55.5084                  0                        3e-10                    6e-10                    8e-10                    2e-10                    9e-10                    -235.731
    S                        5.7e-19                  5.7e-19                  0                        0                        0                        0                        0                        0                        0
    Si                       0.000100001              0.0001                   0                        0                        0                        3e-10                    1e-10                    2e-10                    -390.838
    Z                        0.00017249               0.00017249               0                        0                        0                        0                        0                        0                        -42.4272
    =================================================================================================================================================================================================================================
    Species                  Amount [mol]             Mole Fraction [mol/mol]  Activity Coefficient [-] Activity [-]             Potential [kJ/mol]       Dual Potential [kJ/mol]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Al+++                    9.89996e-08              1.50498e-09              3.93858e-07              3.89923e-14              -557.515                 2.80635e-13
    AlO+                     3.85672e-08              5.86294e-10              0.292194                 1.12693e-08              -708.392                 7.20374e-13
    AlO2-                    9.60796e-05              1.46059e-06              0.292196                 2.80744e-05              -859.268                 2.89164e-16
    AlOH++                   7.36002e-09              1.11886e-10              0.00401324               2.95379e-11              -753.172                 3.77483e-12
    CO(aq)                   1.83254e-10              2.78581e-12              1                        1.83257e-10              -186.558                 1.51608e-10
    CO2(aq)                  3.76868e-06              5.72911e-08              3                        1.13062e-05              -422.289                 7.37202e-15
    CO3--                    1.01827e-06              1.54796e-08              0.0358628                3.65184e-08              -573.166                 2.72844e-14
    Ca(HCO3)+                1.79306e-08              2.72579e-10              0.292194                 5.23929e-09              -1201.58                 1.54946e-12
    Ca(HSiO3)+               1.08558e-10              1.65029e-12              0.292194                 3.17205e-11              -1641.59                 2.55925e-10
    Ca++                     1.9264e-05               2.92849e-07              0.36756                  7.08075e-06              -583.637                 1.44221e-15
    CaCO3(aq)                1.27512e-09              1.93842e-11              1                        1.27514e-09              -1156.8                  2.17884e-11
    CaCl+                    6.36868e-05              9.6816e-07               0.292194                 1.86092e-05              -713.513                 4.36241e-16
    CaCl2(aq)                1.70291e-05              2.58875e-07              1                        1.70293e-05              -843.39                  1.63149e-15
    CaOH+                    8.21547e-10              1.24891e-11              0.292194                 2.40055e-10              -779.294                 3.38176e-11
    CaSO4(aq)                1e-20                    1.52019e-22              1                        1.00001e-20              -1437.8                  0
    Cl-                      4.00233                  0.0608431                0.788966                 3.15775                  -129.876                 6.94164e-21
    ClO-                     1.50591e-22              2.28927e-24              0.292196                 4.40026e-23              -181.116                 184.492
    ClO2-                    5.91892e-23              8.99788e-25              0.292196                 1.72951e-23              -131.95                  469.389
    ClO3-                    4.10539e-23              6.24097e-25              0.292196                 1.19959e-23              -160.33                  676.74
    ClO4-                    3.0521e-23               4.63977e-25              0.292196                 8.91822e-24              -162.516                 910.285
    H+                       5.23784e-08              7.96251e-10              1.90916                  1e-07                    -44.7805                 5.30424e-13
    H2(aq)                   0.00075                  1.14014e-05              1                        0.000750011              -4.70661                 3.70437e-17
    H2O(l)                   55.5077                  0.843822                 0.967394                 0.816309                 -240.438                 5.00521e-22
    H2O2(aq)                 1.40968e-22              2.14298e-24              1                        1.4097e-22               -279.083                 197.086
    H2S(aq)                  1e-20                    1.52019e-22              1                        1.00001e-20              -160.773                 0
    H2S2O3(aq)               1e-20                    1.52019e-22              1                        1.00001e-20              -670.529                 0
    H2S2O4(aq)               1e-20                    1.52019e-22              1                        1.00001e-20              -752.688                 0
    HAlO2(aq)                3.77545e-06              5.7394e-08               1                        3.7755e-06               -904.049                 7.3588e-15
    HCO3-                    9.51723e-05              1.4468e-06               0.523504                 4.98238e-05              -617.946                 2.91921e-16
    HCl(aq)                  5.06501e-08              7.69978e-10              1                        5.06508e-08              -174.657                 5.48524e-13
    HClO(aq)                 1.49649e-22              2.27495e-24              1                        1.49651e-22              -224.735                 185.653
    HClO2(aq)                5.52827e-23              8.40401e-25              1                        5.52834e-23              -143.561                 502.559
    HO2-                     1.26278e-22              1.91967e-24              0.292196                 3.68985e-23              -211.376                 220.013
    HS-                      1e-20                    1.52019e-22              0.292196                 2.922e-21                -121.66                  0
    HS2O3-                   1e-20                    1.52019e-22              0.292196                 2.922e-21                -668.203                 0
    HS2O4-                   1e-20                    1.52019e-22              0.292196                 2.922e-21                -751.617                 0
    HSO3-                    1e-20                    1.52019e-22              0.292196                 2.922e-21                -664.121                 0
    HSO4-                    1e-20                    1.52019e-22              0.522344                 5.22351e-21              -890.08                  0
    HSO5-                    1e-20                    1.52019e-22              0.292196                 2.922e-21                -776.845                 0
    HSiO3-                   1.11635e-06              1.69706e-08              0.292196                 3.26197e-07              -1057.96                 2.48872e-14
    K+                       9.87448e-10              1.50111e-11              0.539955                 5.33185e-10              -345.45                  2.8136e-11
    KCl(aq)                  1.25523e-11              1.90818e-13              1                        1.25524e-11              -475.326                 2.21337e-09
    KHSO4(aq)                1e-20                    1.52019e-22              1                        1.00001e-20              -1155.13                 0
    KOH(aq)                  2.222e-16                3.37787e-18              1                        2.22204e-16              -541.107                 0.000125035
    KSO4-                    1e-20                    1.52019e-22              0.292196                 2.922e-21                -1168.49                 0
    Mg(HCO3)+                2.08489e-08              3.16943e-10              0.292194                 6.09202e-09              -1099.48                 1.33258e-12
    Mg(HSiO3)+               1.85483e-10              2.81969e-12              0.292194                 5.41978e-11              -1539.49                 1.49786e-10
    Mg++                     1.72811e-05              2.62705e-07              0.473801                 8.1879e-06               -481.529                 1.6077e-15
    MgCO3(aq)                4.89943e-10              7.44806e-12              1                        4.89949e-10              -1054.7                  5.67062e-11
    MgCl+                    8.26936e-05              1.2571e-06               0.292194                 2.41629e-05              -611.406                 3.35973e-16
    MgOH+                    3.82816e-09              5.81953e-11              0.591789                 2.26549e-09              -677.186                 7.25748e-12
    MgSO4(aq)                1e-20                    1.52019e-22              1                        1.00001e-20              -1336.91                 0
    Na+                      4.00248                  0.0608453                0.788914                 3.15766                  -260.892                 6.94139e-21
    NaCl(aq)                 2.26749                  0.0344701                1                        2.26752                  -390.768                 1.22527e-20
    NaHSiO3(aq)              3.03392e-05              4.61213e-07              1                        3.03396e-05              -1318.85                 9.15739e-16
    NaOH(aq)                 1.66518e-06              2.53139e-08              1                        1.6652e-06               -456.549                 1.66846e-14
    NaSO4-                   1e-20                    1.52019e-22              0.292196                 2.922e-21                -1145.5                  0
    O2(aq)                   8.11732e-23              1.23399e-24              1                        8.11743e-23              -129.197                 342.265
    OH-                      1.5076e-06               2.29184e-08              0.528428                 7.96672e-07              -195.657                 1.84284e-14
    S2--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -64.2461                 0
    S2O3--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -667.777                 0
    S2O4--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -746.566                 0
    S2O5--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -937.405                 0
    S2O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1113.91                 0
    S2O8--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1266.9                  0
    S3--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -71.5219                 0
    S3O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1106.01                 0
    S4--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -77.5273                 0
    S4O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1192.88                 0
    S5--                     1e-20                    1.52019e-22              0.0040133                4.01335e-23              -82.2767                 0
    S5O6--                   1e-20                    1.52019e-22              0.0040133                4.01335e-23              -1107.1                  0
    SO2(aq)                  1e-20                    1.52019e-22              1                        1.00001e-20              -435.356                 0
    SO3--                    1e-20                    1.52019e-22              0.0040133                4.01335e-23              -628.241                 0
    SO4--                    1e-20                    1.52019e-22              0.0319362                3.19367e-22              -882.113                 0
    SiO2(aq)                 6.85442e-05              1.042e-06                1                        6.85451e-05              -862.3                   4.05327e-16
    Halite                   1e-10                    1                        1                        1                        -386.836                 3.93251
    Calcite                  1e-10                    1                        1                        1                        -1132.71                 24.0922
    Dolomite                 1e-10                    1                        1                        1                        -2172.27                 39.2262
    K-Feldspar               1e-10                    1                        1                        1                        -3754.4                  37.2138
    Quartz                   1e-10                    1                        1                        1                        -857.834                 4.46604
    Kaolinite                1e-10                    1                        1                        1                        -3796.98                 -23.8434
    =================================================================================================================================================================================================================================
    Phase                    Amount [mol]             Stability                Stability Index [-]      Mass [kg]                Volume [m3]              Density [kg/m3]          Molar Volume [m3/mol]    Volume Fraction [m3/m3]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    Aqueous                  65.7813                  stable                   -9.64327e-17             1.36645                  0.00117003               1167.87                  1.77867e-05              1
    Halite                   1e-10                    under stable             -0.614721                5.84428e-12              2.7015e-15               2163.35                  2.7015e-05               2.30891e-12
    Calcite                  1e-10                    under stable             -3.76604                 1.00087e-11              3.6934e-15               2709.89                  3.6934e-05               3.15666e-12
    Dolomite                 1e-10                    under stable             -6.13176                 1.84401e-11              6.4365e-15               2864.92                  6.4365e-05               5.50112e-12
    K-Feldspar               1e-10                    under stable             -5.81718                 2.78332e-11              1.0887e-14               2556.55                  0.00010887               9.30486e-12
    Quartz                   1e-10                    under stable             -0.698121                6.00843e-12              2.2688e-15               2648.29                  2.2688e-05               1.93909e-12
    Kaolinite                1e-10                    super stable              3.72715                 2.5816e-11               9.952e-15                2594.06                  9.952e-05                8.50574e-12
    =================================================================================================================================================================================================================================
    Ionic Strength [molal]   pH                       pE                       Reduction Potential [V]  Alkalinity [eq/L]
    ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    4.00271                  7                        -0.284906                -0.01889                 0.000159835
    =================================================================================================================================================================================================================================
    */

    // Differences         of    PHREEQC   vs Reaktoro
    // Ionic Strength [molal]    6.27         4.00271
    // Total alkalinity (eq/kg)  4.949e-04    0.000159835
    // # iterations              11           74
    // pe                        4.0          -0.284906
    // Because of the added water in Reaktoro:
    // H(0)                      2.711e-25    111.017
    // O(0)                      5.309e-34    55.5084
    // even though
    // H2O                       5.553e+01    55.5077
    //
    // Which species can be crucial for the comparison?

    state0.setSpeciesAmount("Calcite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Quartz", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Halite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Dolomite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("Kaolinite", initial_mineral_amount_mol, "mol");
    state0.setSpeciesAmount("K-Feldspar", initial_mineral_amount_mol, "mol");

    // Set initial value of mineral states for all the reactions
//    for (Reaction& reaction : reaction_minerals)
//        reaction.setInitialAmounts(state0.speciesAmounts());
//    for(auto i = 0; i < reaction_minerals.size(); i++)
//        reaction_minerals[i].setInitialAmounts(state0.speciesAmounts());
    reaction_halite.setInitialAmounts(state0.speciesAmounts());
    reaction_calcite.setInitialAmounts(state0.speciesAmounts());
    reaction_dolomite.setInitialAmounts(state0.speciesAmounts());
    reaction_quartz.setInitialAmounts(state0.speciesAmounts());
    reaction_kfeldspar.setInitialAmounts(state0.speciesAmounts());
    reaction_kaolinite.setInitialAmounts(state0.speciesAmounts());

    KineticPath path(reactions, partition);

    ChemicalOutput output = path.output();
    output.filename(result_filename);
    output.add("time(units=second)");
    output.add("speciesMolality(Na+ units=mmolal)", "Na+ [mmol]");
    output.add("speciesMolality(Cl- units=mmolal)", "Cl- [mmol]");
    output.add("speciesMolality(Mg++ units=mmolal)", "Mg++ [mmol]");
    output.add("speciesMolality(Ca++ units=mmolal)", "Ca++ [mmol]");
    output.add("speciesMolality(CO3-- units=mmolal)", "CO3-- [mmol]");
    output.add("speciesMolality(K+ units=mmolal)", "K+ [mmol]");
    output.add("speciesMolality(Al+++ units=mmolal)", "Al+++ [mmol]");
    output.add("speciesMolality(Calcite units=molal)", "Calcite [mol]");
    output.add("speciesMolality(Dolomite units=molal)", "Dolomite [mol]");
    output.add("speciesMolality(Halite units=molal)", "Halite [mol]");
    output.add("speciesMolality(K-Feldspar units=molal)", "K-Feldspar [mol]");
    output.add("speciesMolality(Quartz units=molal)", "Quartz [mol]");
    output.add("speciesMolality(Kaolinite units=molal)", "Kaolinite [mol]");
    output.add("pH");

    tic(KINETIC_PATH);

    ChemicalState state0_tmp = state0;

    path.solve(state0, t0, dt, n, "second");

    auto kinetic_path_time = toc(KINETIC_PATH);
    std::cout << "total kinetic path time = " << kinetic_path_time << " s" << std::endl;

    //std::cout << "state0 = \n" << state0 << std::endl;

    double delta_halite = state0.speciesAmount("Halite") - initial_mineral_amount_mol;
    double delta_calcite = state0.speciesAmount("Calcite") - initial_mineral_amount_mol;
    double delta_dolomite = state0.speciesAmount("Dolomite") - initial_mineral_amount_mol;
    double delta_kfeldspar = state0.speciesAmount("K-Feldspar") - initial_mineral_amount_mol;
    double delta_quartz = state0.speciesAmount("Quartz") - initial_mineral_amount_mol;
    double delta_kaolinite = state0.speciesAmount("Kaolinite") - initial_mineral_amount_mol;

    std::cout << "delta_halite = " << delta_halite << std::endl;
    std::cout << "delta_calcite = " << delta_calcite << std::endl;
    std::cout << "delta_dolomite = " << delta_dolomite << std::endl;
    std::cout << "delta_kfeldspar = " << delta_kfeldspar << std::endl;
    std::cout << "delta_quartz = " << delta_quartz << std::endl;
    std::cout << "delta_kaolinite = " << delta_kaolinite << std::endl;

    //    Differences       PHREEQC   vs Reaktoro
    // --------------------------------------------------
    //    Halite           -0.2442      -0.624908
    //    Calcite          -4.359e-06   -4.38649e-06
    //    Dolomite-dis     -1.068e-06   -1.22678e-06
    //    K-Feldspar       -2.139e-12   -2.14051e-12
    //    Quartz           -1.776e-12   -2.14406e-12
    //    Kaolinite         1.361e-09    7.10543e-15
}
