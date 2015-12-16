//
// Created by Matthew Guay on 9/17/15.
//

#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cassert>

#include "cell_cpp.h"

namespace kcnet {

    const double V_offset = 65.;

    CppKenyonCell::CppKenyonCell(double i_dt) {
        dt = i_dt; // RK4 update time-step

        V = -69.145387; // membrane voltage (mV)
        C = 2.9e-1; // membrane capacitance (μF) times V-> mV conversion factor (1000)
        Cinv = 1. / C; // reciprocal of membrane capacitance
        Ca = 2.629407e-4; // intracellular Ca2+ concentration

        // Leakage current parameters.
        g_L = 2.9e-2; // leakage conductance (μS)
        g_KL = 1.16e-3; // potassium leakage conductance (μS)
        E_L = -65.; // leakage reversal potential (mV)

        // Transient calcium current parameters (I_Ca).
        g_Ca = 2.9e-2; // I_Ca conductance (μS)
        E_Ca = 114.390135; // I_Ca reversal potential (mV)
        m_Ca = 8.972934e-2; // I_Ca m gating variable
        h_Ca = 1.023953e-1; // I_Ca h gating variable

        // Ca-dependent potassium current (I_KCa).
        g_KCa = 2.9e-1; // I_KCa conductance (μS)
        E_KCa = -95.; // I_CKa reversal potential (mV)
        m_KCa = 2.304748e-4; // I_KCa m gating variable

        // Transient (A-type) potassium current (I_KA).
        g_KA = 1.45e-2; // I_KA conductance (μS)
        E_KA = -95.; // I_KA reversal potential (mV)
        m_KA = 2.542774e-1; // I_KA m gating variable

        // Sodium current (I_Na).
        g_Na = 26.1; // I_Na conductance (μS)
        E_Na = 50.; // I_Na reversal potential (mV)
        m_Na = 6.151612e-3; // I_Na m gating variable
        h_Na = 9.985889e-1; // I_Na h gating variable

        // Potassium current (I_K).
        g_K = 2.9; // I_K conductance (μS)
        E_K = -95.; // I_K reversal potential (mV)
        m_K = 1.875148e-2; // I_K m gating variable
        h_K = 9.546243e-1; // I_K h gating variable

        // Current info.
        I_L = -1.139e-1;
        I_KL = 3.025e-2;
        I_Ca = -1.1435e-2;
        I_KCa = 8.226e-7;
        I_KA = 9.799e-2;
        I_Na = -2.940e-8;
        I_K = 1.513e-10;
        I_syn = 0;

        weight[0] = 0.;
        weight[1] = 0.5 * dt;
        weight[2] = 0.5 * dt;
        weight[3] = dt;
    }

    CppKenyonCell::~CppKenyonCell(){}

    void CppKenyonCell::update(int i_n_synapses, std::vector<double> i_syn_gOcs, double i_syn_E, double i_I_noise) {
        /*double k_V[4], k_Ca[4], k_m_Ca[4], k_h_Ca[4], k_m_KCa[4],
               k_m_KA[4], k_m_Na[4], k_h_Na[4], k_m_K[4], k_h_K[4];*/
        double k_V[4] = {0., 0., 0., 0.};
        double k_Ca[4] = {0., 0., 0., 0.};
        double k_m_Ca[4] = {0., 0., 0., 0.};
        double k_h_Ca[4] = {0., 0., 0., 0.};
        double k_m_KCa[4] = {0., 0., 0., 0.};
        double k_m_KA[4] = {0., 0., 0., 0.};
        double k_m_Na[4] = {0., 0., 0., 0.};
        double k_h_Na[4] = {0., 0., 0., 0.};
        double k_m_K[4] = {0., 0., 0., 0.};
        double k_h_K[4] = {0., 0., 0., 0.};

        // Each of the four iterations is a step in the RK4 update procedure.
        double I_syn_i = 0; // Declare at higher scope, it's useful to use its last value to record I_syn.
        for (int i = 0; i < 4; i++) {
            double w_i, V_i, Ca_i, m_Ca_i, h_Ca_i, m_KCa_i, m_KA_i,
                    m_Na_i, h_Na_i, m_K_i, h_K_i;
            double I_L_i, I_KL_i, I_Ca_i, I_KCa_i, I_KA_i,
                    I_Na_i, I_K_i;

            // Get circular index of previous step's i value.
            int iprev;
            (i == 0) ? (iprev = 3) : (iprev = i - 1);

            // This RK4 step's weight.
            w_i = weight[i];

            // Update this RK4 step's variable values.
            V_i = V + w_i * k_V[iprev];
            Ca_i = Ca + w_i * k_Ca[iprev];
            m_Ca_i = m_Ca + w_i * k_m_Ca[iprev];
            h_Ca_i = h_Ca + w_i * k_h_Ca[iprev];
            m_KCa_i = m_KCa + w_i * k_m_KCa[iprev];
            m_KA_i = m_KA + w_i * k_m_KA[iprev];
            m_Na_i = m_Na + w_i * k_m_Na[iprev];
            h_Na_i = h_Na + w_i * k_h_Na[iprev];
            m_K_i = m_K + w_i * k_m_K[iprev];
            h_K_i = h_K + w_i * k_h_K[iprev];

            // Update this RK4 step's I_Ca reversal potential.
            E_Ca = 12.8 * std::log(2. / Ca_i);

            // This RK4 step's current calculations.
            I_L_i = g_L * (V_i - E_L);
            I_KL_i = g_KL * (V_i - E_K);
            I_Ca_i = g_Ca * std::pow(m_Ca_i, 2) * h_Ca_i * (V_i - E_Ca);
            I_KCa_i = g_KCa * std::pow(m_KCa_i, 2) * (V_i - E_KCa);
            I_KA_i = g_KA * m_KA_i * (V_i - E_KA);
            I_Na_i = g_Na * std::pow(m_Na_i, 3) * h_Na_i * (V_i - E_Na);
            I_K_i = g_K * std::pow(m_K_i, 4) * h_K_i * (V_i - E_K);

            // Calculate synaptic currents.
            I_syn_i = 0;
            if (i_n_synapses > 0) {
                for (int j = 0; j < i_n_synapses; j++) {
                    I_syn_i += i_syn_gOcs[i * i_n_synapses + j] * (V_i - i_syn_E);
                }
            }

            // Update RK4 k variables.
            k_V[i] = Cinv * (-I_L_i - I_KL_i - I_Ca_i - I_KCa_i - I_KA_i - I_Na_i - I_K_i - I_syn_i - i_I_noise);
            k_Ca[i] = -5.2e-5 * I_Ca_i - (Ca_i - 2.4e-4) / 100.;
            k_m_Ca[i] = (minf_Ca(V_i) - m_Ca_i) / mtau_Ca(V_i);
            k_h_Ca[i] = (hinf_Ca(V_i) - h_Ca_i) / htau_Ca(V_i);
            k_m_KCa[i] = (minf_KCa(Ca_i) - m_KCa_i) / mtau_KCa(Ca_i);
            k_m_KA[i] = (minf_KA(V_i) - m_KA_i) / mtau_KA(V_i);
            k_m_Na[i] = malpha_Na(V_i) * (1. - m_Na_i) - mbeta_Na(V_i) * m_Na_i;
            k_h_Na[i] = halpha_Na(V_i) * (1. - h_Na_i) - hbeta_Na(V_i) * h_Na_i;
            k_m_K[i] = malpha_K(V_i) * (1. - m_K_i) - mbeta_K(V_i) * m_K_i;
            k_h_K[i] = halpha_K(V_i) * (1. - h_K_i) - hbeta_K(V_i) * h_K_i;
        }

        // Update state variables.
        V += dt/6 * (k_V[0] + 2*k_V[1] + 2*k_V[2] + k_V[3]);
        Ca += dt/6 * (k_Ca[0] + 2*k_Ca[1] + 2*k_Ca[2] + k_Ca[3]);
        m_Ca += dt/6 * (k_m_Ca[0] + 2*k_m_Ca[1] + 2*k_m_Ca[2] + k_m_Ca[3]);
        h_Ca += dt/6 * (k_h_Ca[0] + 2*k_h_Ca[1] + 2*k_h_Ca[2] + k_h_Ca[3]);
        m_KCa += dt/6 * (k_m_KCa[0] + 2*k_m_KCa[1] + 2*k_m_KCa[2] + k_m_KCa[3]);
        m_KA += dt/6 * (k_m_KA[0] + 2*k_m_KA[1] + 2*k_m_KA[2] + k_m_KA[3]);
        m_Na += dt/6 * (k_m_Na[0] + 2*k_m_Na[1] + 2*k_m_Na[2] + k_m_Na[3]);
        h_Na += dt/6 * (k_h_Na[0] + 2*k_h_Na[1] + 2*k_h_Na[2] + k_h_Na[3]);
        m_K += dt/6 * (k_m_K[0] + 2*k_m_K[1] + 2*k_m_K[2] + k_m_K[3]);
        h_K += dt/6 * (k_h_K[0] + 2*k_h_K[1] + 2*k_h_K[2] + k_h_K[3]);

        // Record current values.
        I_L = g_L * (V - E_L);
        I_KL = g_KL * (V - E_K);
        I_Ca = g_Ca * std::pow(m_Ca, 2) * h_Ca * (V - E_Ca);
        I_KCa = g_KCa * std::pow(m_KCa, 2) * (V - E_KCa);
        I_KA = g_KA * m_KA * (V - E_KA);
        I_Na = g_Na * std::pow(m_Na, 3) * h_Na * (V - E_Na);
        I_K = g_K * std::pow(m_K, 4) * h_K * (V - E_K);
        I_syn = I_syn_i;

        // Initialise the reset_to_this state to the initial state here.
        reset_to_this.push_back(V);
        reset_to_this.push_back(Ca);
        reset_to_this.push_back(m_Ca);
        reset_to_this.push_back(h_Ca);
        reset_to_this.push_back(m_KCa);
        reset_to_this.push_back(m_KA);
        reset_to_this.push_back(m_Na);
        reset_to_this.push_back(h_Na);
        reset_to_this.push_back(m_K);
        reset_to_this.push_back(h_K);
        reset_to_this.push_back(I_L);
        reset_to_this.push_back(I_KL);
        reset_to_this.push_back(I_Ca);
        reset_to_this.push_back(I_KCa);
        reset_to_this.push_back(I_KA);
        reset_to_this.push_back(I_Na);
        reset_to_this.push_back(I_K);
        reset_to_this.push_back(I_syn);
    }

    // I_Ca gating variable functions (taken from Maxim Bazhenov's RE cell code).
    double CppKenyonCell::minf_Ca(double i_V) {
        return 1. / (1. + std::exp(-(i_V + 52.) / 7.4));
    }
    double CppKenyonCell::mtau_Ca(double i_V) {
        return (3. + 1. / (std::exp((i_V + 27.) / 10.) + std::exp(-(i_V + 102.) / 15.))) / 9.90;
    }
    double CppKenyonCell::hinf_Ca(double i_V) {
        return 1. / (1. + std::exp((i_V + 80.) / 5.));
    }
    double CppKenyonCell::htau_Ca(double i_V) {
        return (85. + 1. / (std::exp((i_V + 48.) / 4.) + std::exp(-(i_V + 407.) / 50.))) / 3.74;
    }

    // I_KCa gating variable functions.
    double CppKenyonCell::minf_KCa(double i_Ca) {
        return 3333. * std::pow(i_Ca, 2) / (3333. * std::pow(i_Ca, 2) + 1.);
    }
    double CppKenyonCell::mtau_KCa(double i_Ca) {
        return std::max(0.1, 0.896 / (100. * std::pow(i_Ca, 2) + 0.03));
    }

    // I_KA gating variable functions.
    double CppKenyonCell::minf_KA(double i_V) {
        return 1. / (1. + std::exp(-(i_V + 60.) / 8.5));
    }
    double CppKenyonCell::mtau_KA(double i_V) {
        return (1. / (std::exp((i_V + 35.82) / 19.69) + std::exp(-(i_V + 79.69) / 12.7)) + 0.37) / 3.74;
    }

    // I_Na gating variable functions.
    double CppKenyonCell::malpha_Na(double i_V) {
        double V2 = i_V + V_offset;
        return std::abs(V2 - 13.) < 1e-9? 1.28 : 0.32 * (13. - V2) / (std::exp((13. - V2) / 4.) - 1.);
    }
    double CppKenyonCell::mbeta_Na(double i_V) {
        double V2 = i_V + V_offset;
        return std::abs(V2 - 40.) < 1e-9? 1.4 : 0.28 * (V2 - 40.) / (std::exp((V2 - 40.) / 5.) - 1.);
    }
    double CppKenyonCell::halpha_Na(double i_V) {
        double V2 = i_V + V_offset;
        return 0.128 * std::exp((17. - V2) / 18.);
    }
    double CppKenyonCell::hbeta_Na(double i_V) {
        double V2 = i_V + V_offset;
        return 4. / (std::exp((40. - V2) / 5.) + 1.);
    }

    // I_K gating variable functions.
    double CppKenyonCell::malpha_K(double i_V) {
        double V2 = i_V + V_offset;
        return std::abs(V2 - 15.) < 1e-9? 0.16 : 0.032 * (15. - V2) / (std::exp((15. - V2) / 5.) - 1.);
    }
    double CppKenyonCell::mbeta_K(double i_V) {
        double V2 = i_V + V_offset;
        return 0.5 * std::exp((10. - V2) / 40.);
    }
    double CppKenyonCell::halpha_K(double i_V) {
        double V2 = i_V + V_offset;
        return 0.028 * std::exp((15. - V2) / 15.) + 2. / (std::exp((85. - V2) / 10.) + 1.);
    }
    double CppKenyonCell::hbeta_K(double i_V) {
        double V2 = i_V + V_offset;
        return 0.4 / (std::exp((40. - V2) / 10.) + 1.);
    }

    void CppKenyonCell::save_state(std::string i_name) {
        try {
            std::ofstream file;
            std::string file_name = "../states/" + i_name + ".txt";
            file.open(file_name);

            file << V << ","
            << Ca << ","
            << m_Ca << ","
            << h_Ca << ","
            << m_KCa << ","
            << m_KA << ","
            << m_Na << ","
            << h_Na << ","
            << m_K << ","
            << h_K << ","
            << I_L << ","
            << I_KL << ","
            << I_Ca << ","
            << I_KCa << ","
            << I_KA << ","
            << I_Na << ","
            << I_K << ","
            << I_syn << "\n";

            file.close();
        } catch(int e) {
            std::cout << "CppKenyonCell state save failed. Exception " << e << "\n";
        }
    }

    void CppKenyonCell::load_state(std::string _name) {

        // Get the location of the file. Open a read stream.
        std::string file_name = "../states/" + _name + ".dat";
        std::ifstream file(file_name);

        if(file.is_open()) {
            // If the file opened successfully:

            // Erase the previous reset state.
            reset_to_this.erase(reset_to_this.begin(), reset_to_this.end());

            // Get each line of data, convert to double and add to vector reset_to_this.
            std::string data;
            while(std::getline(file, data, ',')) {
                reset_to_this.push_back(std::stod(data));
            }
            file.close();
        } else {
            // File could not be opened.
            printf("Warning: could not load state %s.", _name.c_str());
        }
    }

    void CppKenyonCell::reset_state(std::vector<int> _idxs) {

        // Tedious but what is there to do
        if(_idxs[0] != 0) { // Load new state for V.
            V = reset_to_this[0];
        }
        if(_idxs[1] != 0) { // Load new state for Ca.
            Ca = reset_to_this[1];
        }
        if(_idxs[2] != 0) { // Load new state for m_Ca.
            m_Ca = reset_to_this[2];
        }
        if(_idxs[3] != 0) { // Load new state for h_Ca.
            h_Ca = reset_to_this[3];
        }
        if(_idxs[4] != 0) { // Load new state for m_KCa.
            m_KCa = reset_to_this[4];
        }
        if(_idxs[5] != 0) { // Load new state for m_KA.
            m_KA = reset_to_this[5];
        }
        if(_idxs[6] != 0) { // Load new state for m_Na.
            m_Na = reset_to_this[6];
        }
        if(_idxs[7] != 0) { // Load new state for h_Na.
            h_Na = reset_to_this[7];
        }
        if(_idxs[8] != 0) {  // Load new state for m_K.
            m_K = reset_to_this[8];
        }
        if(_idxs[9] != 0) { // Load new state for h_K.
            h_K = reset_to_this[9];
        }
        if(_idxs[10] != 0) { // Load new state for I_L.
            I_L = reset_to_this[10];
        }
        if(_idxs[11] != 0) { // Load new state for I_KL.
            I_KL = reset_to_this[11];
        }
        if(_idxs[12] != 0) { // Load new state for I_Ca.
            I_Ca = reset_to_this[12];
        }
        if(_idxs[13] != 0) {
            I_KCa = reset_to_this[13];
        }
        if(_idxs[14] != 0) {
            I_KA = reset_to_this[14];
        }
        if(_idxs[15] != 0) {
            I_Na = reset_to_this[15];
        }
        if(_idxs[16] != 0) {
            I_K = reset_to_this[16];
        }
        if(_idxs[17] != 0) {
            I_syn = reset_to_this[17];
        }
    }

};

