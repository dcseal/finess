#ifndef MAXWELLPARAMS_H
#define MAXWELLPARAMS_H

// TODO - move this into global_ini_params - right now there is no call that
// will set this up in the code!!! -DS

struct MaxwellParams{

 private:

   // light speed
   double cs_light;
   double cs_light_squared;

   // correction potential speed
   double cp_speed;
   double cp_speed_squared;

   // correction potential decay rates
   double eps_Ecp; // decay rate of electric field correction potential
   double eps_Bcp; // decay rate of magnetic field correction potential
   int clean_E_field;

 public:

   const double get_cs_light()const{return cs_light;}
   const double get_cs_light_squared()const{return cs_light_squared;}
   const double get_one_over_epsilon()const{return cs_light_squared;}
   const double get_cp_speed()const{return cp_speed;}
   const double get_cp_speed_squared()const{return cp_speed_squared;}
   const double get_eps_Ecp()const{return eps_Ecp; }
   const double get_eps_Bcp()const{return eps_Bcp; }
   const int get_clean_E_field()const{
     return clean_E_field; }

 public:
   MaxwellParams():
      cs_light(0),
      cs_light_squared(0),
      cp_speed(0),
      cp_speed_squared(0),
      eps_Ecp(0),
      eps_Bcp(0),
      clean_E_field(0){}
   void set_cs_light(double in){
     cs_light = in;
     cs_light_squared = in*in;}
   void set_cp_speed(double in){
     cp_speed = in;
     cp_speed_squared = in*in; }
   void set_eps_Ecp(double in){eps_Ecp = in;}
   void set_eps_Bcp(double in){eps_Bcp = in;}
   void set_clean_E_field(int in){clean_E_field = in;}

};
extern MaxwellParams maxwellParams;

#endif
