#ifndef _STATEVARS_H_
#define _STATEVARS_H_

#include <algorithm>

#include "tensors.h"


class StateVars{
    private:
        double t;
        dTensorBC2 q;
        dTensorBC2 aux;
    public:
        StateVars(double t, int mx, int meqn, int maux, int mbc ):            
            t(t),
            q(mx, meqn, mbc),
            aux(mx, std::max(1, maux), mbc)
        { }

        double get_t() const{
            return this->t;
        }

        void set_t(double t){
            this->t = t;
        }
        
        const dTensorBC2& const_ref_q() const{
            return this->q;
        }

        dTensorBC2& ref_q(){
            return this->q;
        }

        const dTensorBC2& const_ref_aux() const{
            return this->aux;
        }

        dTensorBC2& ref_aux(){
            return this->aux;
        }

};

#endif
