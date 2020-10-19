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

#include <algorithm>

#include <Reaktoro/Reaktoro.hpp>
using namespace Reaktoro;

int main()
{
    double t0 = 0;
    double dt = 0.1; // in minutes
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

    //std::string activity_model = "pitzer";
    std::string activity_model = "hkf";
    std::cout << "Activity model: " << activity_model << std::endl;


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

    std::vector<std::string> minerals = {"Halite", "Calcite", "Dolomite", "K-Feldspar", "Kaolinite", "Quartz"};
    auto num_minerals = minerals.size();

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
    std::string eq_str_quartz = "Quartz + H2O(l) = H+ + HSiO3-";
    //std::string eq_str_quartz = "Quartz = SiO2(aq)";
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
    ReactionSystem reactions(system, {reaction_halite, reaction_calcite, reaction_dolomite, reaction_kfeldspar, reaction_kaolinite, reaction_quartz});

    std::vector<Partition> partitions(num_minerals + 1);
    std::fill(partitions.begin(), partitions.end(), Partition(system));

    for(int i = 0; i < num_minerals + 1; i++){
        auto kinetic_minerals = std::vector<std::string>(minerals.begin() + i, minerals.end());
        std::cout << "Kinetic minerals: "; for(const auto& mineral : kinetic_minerals) std::cout << mineral << ", "; std::cout << std::endl;

        // Initialize current partition with kinetic minerals
        partitions[i].setKineticPhases(kinetic_minerals);
        auto num_kinetic_minerals = kinetic_minerals.size();

        double T = 61.0 + 273.15;       // temperature (in units of celsius)
        double P = 1 * 1.01325 * 1e5;   // pressure (in units of atm)
        EquilibriumInverseProblem problem(partitions[i]);
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


        EquilibriumOptions options;
        options.hessian = GibbsHessian::ApproximationDiagonal;

        EquilibriumInverseSolver solver(partitions[i]);
        solver.setOptions(options);
        ChemicalState state0(system);
        EquilibriumResult result = solver.solve(state0, problem);

        std::cout << "# iter (in initial state) = " << result.optimum.iterations << std::endl;

        double initial_mineral_amount_mol = 10.0;
        state0.setSpeciesAmount("Calcite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Quartz", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Halite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Dolomite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("Kaolinite", initial_mineral_amount_mol, "mol");
        state0.setSpeciesAmount("K-Feldspar", initial_mineral_amount_mol, "mol");

        reaction_halite.setInitialAmounts(state0.speciesAmounts());
        reaction_calcite.setInitialAmounts(state0.speciesAmounts());
        reaction_dolomite.setInitialAmounts(state0.speciesAmounts());
        reaction_kfeldspar.setInitialAmounts(state0.speciesAmounts());
        reaction_kaolinite.setInitialAmounts(state0.speciesAmounts());
        reaction_quartz.setInitialAmounts(state0.speciesAmounts());

        KineticOptions kinetic_options;
        kinetic_options.equilibrium = options;

        KineticPath path(reactions, partitions[i]);
        path.setOptions(kinetic_options);

        std::string result_filename = "kinetics-benchmark-kin-" + std::to_string(num_kinetic_minerals)
                                                           + "-eq-" + std::to_string(num_minerals - num_kinetic_minerals)
                                                           + "-dt-" + std::to_string(dt)
                                                           + "-n-" +  std::to_string(n) +
                                                           + "-ncells-" +  std::to_string(ncells) +
                                                           + "-activity-model-" + activity_model +
                                                           + "-uniform-sa-1.0-m2mol.txt";
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

        path.solve(state0, t0, dt, n, "second");

        auto kinetic_path_time = toc(KINETIC_PATH);
        std::cout << "----------------------------------------" << std::endl;
        std::cout << "# kin : " << num_kinetic_minerals
                  << ", # eq  : " << num_minerals - num_kinetic_minerals
                  << ", time  : " << kinetic_path_time << " s" << std::endl;
        std::cout << "----------------------------------------" << std::endl;

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

        std::cout << "*******************************************************************************" << std::endl;
    }

}
