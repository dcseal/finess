#ifndef _STATEVARS_H_
#define _STATEVARS_H_

#include <algorithm>

#include "tensors.h"


class StateVars{
    private:
        double t;
        dTensorBC3 q;
        dTensorBC3 aux;
    public:
        StateVars(double t, int mx, int my, int meqn, int maux, int mbc ):            
            t(t),
            q(mx, my, meqn, mbc),
            aux(mx, my, std::max(1, maux), mbc)
        { }

        double get_t() const{
            return this->t;
        }

        void set_t(double t){
            this->t = t;
        }
        
        const dTensorBC3& const_ref_q() const{
            return this->q;
        }

        dTensorBC3& ref_q(){
            return this->q;
        }

        const dTensorBC3& const_ref_aux() const{
            return this->aux;
        }

        dTensorBC3& ref_aux(){
            return this->aux;
        }

};

#endif
