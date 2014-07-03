// Including this header file allows access to all the 1D paramters described in
// the [grid] section of any 1D parameters.ini file.
//
// The singleton dogParamsCart1 allows access to these parameters.  For example,
// after including DogParamsCart1.h, one may accesses mx, the number of grid cells
// via:
//
//      dogParamsCart3.get_mx()
//
// See also: DogParams.h, DogParamsCart3.h, IniDocument.h (global parser)

#ifndef DOGPARAMSCART1_H
#define DOGPARAMSCART1_H

class IniDocument;   // The parser used in DoGPack for reading parameters.ini

struct DogParamsCart3
{

    // -- Methods -- //
    public:

        bool get_is_initialized(){return is_initialized;}
        DogParamsCart3(){is_initialized=false;}

        void init(IniDocument& ini_doc);
        void append_qhelp(const char* filename);

        // User feedback for the welcome screen
        void reportParameters();

        // Accessors
        //
        // Note: melems == mx for backwards compatability
        const int   & get_mx()      const{ return  mx;   }
        const int   & get_my()      const{ return  my;   }
        const int   & get_mz()      const{ return  mz;   }

        const int   & get_mbc()     const{ return  mbc;  }

        const double& get_xlow()    const{ return  xlow; }
        const double& get_xhigh()   const{ return  xhigh;}
        const double& get_dx()      const{ return  dx;   }

        const double& get_ylow()    const{ return  ylow; }
        const double& get_yhigh()   const{ return  yhigh;}
        const double& get_dy()      const{ return  dy;   }

        const double& get_zlow()    const{ return  zlow; }
        const double& get_zhigh()   const{ return  zhigh;}
        const double& get_dz()      const{ return  dz;   }

        // Setters
        void set_mx(int mx_in){ mx = mx_in; }
        void set_my(int my_in){ my = my_in; }
        void set_mz(int mz_in){ mz = mz_in; }

        void set_xlims(double,double);
        void set_ylims(double,double);
        void set_zlims(double,double);

    private:
        void checkParameters();
        void setDerivedParameters();

    // -- Fields -- //
    private:

        bool is_initialized;    // flag defining whether or not the section has been read
        int mx;                 // Number of elements in x-direction
        int my;                 // Number of elements in y-direction
        int mz;                 // Number of elements in z-direction
        int mbc;                // number of ghost cells on each side

        double xlow;            // Left end of domain
        double xhigh;           // Right end of domain
        double dx;              // Derived from other parameters

        double ylow;            // Bottom end of domain
        double yhigh;           // Top end of domain
        double dy;              // Derived from other parameters

        double zlow;            // Bottom end of domain
        double zhigh;           // Top end of domain
        double dz;              // Derived from other parameters

};

extern DogParamsCart3 dogParamsCart3;

#endif
