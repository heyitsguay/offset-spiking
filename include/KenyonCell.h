//
// Created by Matthew Guay on 9/17/15.
//

#ifndef OFFSET_SPIKING_KENYON_CELL_H
#define OFFSET_SPIKING_KENYON_CELL_H
#pragma once

#include <vector>

namespace kcnet{

    const std::vector<int> idxs_all = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    const std::vector<int> idxs_noI = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0};

    class KenyonCell{
    private:
        double weight [4];
        std::vector<double> reset;
    public:
        double dt, V, C, Cinv, Ca,
                g_L, g_KL, E_L, g_Ca, E_Ca,
                m_Ca, h_Ca, g_KCa, E_KCa, m_KCa,
                g_KA, E_KA, m_KA, g_Na, E_Na,
                m_Na, h_Na, g_K, E_K, m_K, h_K,
                I_L, I_KL, I_Ca, I_KCa, I_KA, I_Na, I_K, I_syn;

        KenyonCell(double dt_);
        ~KenyonCell();

        void update(int n_synapses_, std::vector<double> syn_gOcs_, double syn_E_, double I_noise_, double I_in_=0.);

        double minf_Ca(double V_);
        double mtau_Ca(double V_);
        double hinf_Ca(double V_);
        double htau_Ca(double V_);

        double minf_KCa(double Ca_);
        double mtau_KCa(double Ca_);

        double minf_KA(double V_);
        double mtau_KA(double V_);

        double malpha_Na(double V_);
        double mbeta_Na(double V_);
        double halpha_Na(double V_);
        double hbeta_Na(double V_);

        double malpha_K(double V_);
        double mbeta_K(double V_);
        double halpha_K(double V_);
        double hbeta_K(double V_);

        void save_state(std::string name_);
        void load_state(std::string name_);
        void reset_state(std::vector<int> idxs_=idxs_all);
    };
}

#endif // OFFSET_SPIKING_KENYON_CELL_H
	    
