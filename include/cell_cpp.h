//
// Created by Matthew Guay on 9/17/15.
//

#pragma once

#include <vector>

namespace kcnet{

    const std::vector<int> idxs_all = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    const std::vector<int> idxs_noI = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};

    class CppKenyonCell{
    private:
        double weight [4];
        std::vector<double> reset_to_this;
    public:
        double dt, V, C, Cinv, Ca,
                g_L, g_KL, E_L, g_Ca, E_Ca,
                m_Ca, h_Ca, g_KCa, E_KCa, m_KCa,
                g_KA, E_KA, m_KA, g_Na, E_Na,
                m_Na, h_Na, g_K, E_K, m_K, h_K,
                I_L, I_KL, I_Ca, I_KCa, I_KA, I_Na, I_K, I_syn;

        CppKenyonCell(double i_dt);
        ~CppKenyonCell();

        void update(int i_n_synapses, std::vector<double> i_syn_gOcs, double i_syn_E, double i_I_noise);

        double minf_Ca(double i_V);
        double mtau_Ca(double i_V);
        double hinf_Ca(double i_V);
        double htau_Ca(double i_V);

        double minf_KCa(double i_Ca);
        double mtau_KCa(double i_Ca);

        double minf_KA(double i_V);
        double mtau_KA(double i_V);

        double malpha_Na(double i_V);
        double mbeta_Na(double i_V);
        double halpha_Na(double i_V);
        double hbeta_Na(double i_V);

        double malpha_K(double i_V);
        double mbeta_K(double i_V);
        double halpha_K(double i_V);
        double hbeta_K(double i_V);

        void save_state(std::string i_name);
        void load_state(std::string i_name);
        void reset_state(std::vector<int> i_idxs=idxs_all);
    };
}
	    
