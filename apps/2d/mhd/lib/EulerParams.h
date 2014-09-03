#ifndef _EULERPARAMS_H_
#define _EULERPARAMS_H_

class IniDocument;
class EulerParams
{   

    public:
        double gamma;   // gass constant. Used for every application
        double x0;      // used for double-mach reflection problem
        int    OPT;     // used for Riemann test problem
        void init(IniDocument& ini_doc);
        void write_eulerhelp(const char* filename);

};
extern EulerParams eulerParams;

#endif
