#ifndef _STATEVARS_H_
#define _STATEVARS_H_

#include <algorithm>
#include "tensors.h"

// This class holds the necessary information for defining all of the state
// variables: ( t, q, aux ) for the 1D case.
class StateVars{
    public:
        StateVars(double t, int mx, int meqn, int maux, int mbc ):            
            t(t),
            q(mx, meqn, mbc),
            aux(mx, maux, mbc)
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

        // Copy the contents from another state variable to this state
        // variable
        void copyfrom( const StateVars& Qin )
        {
            this->q.copyfrom( Qin.const_ref_q() );
            this->aux.copyfrom( Qin.const_ref_aux() );
            this->t = Qin.get_t();
        }

    private:
        double t;
        dTensorBC2 q;
        dTensorBC2 aux;

};

#endif
