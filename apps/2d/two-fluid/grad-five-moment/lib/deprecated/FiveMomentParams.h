#ifndef FIVEMOMENTPARAMS_H
#define FIVEMOMENTPARAMS_H

struct FiveMomentParams{
   double gamma;
 public:
   const double& get_gamma(){return gamma;}
};

extern FiveMomentParams fiveMomentParams;

#endif
