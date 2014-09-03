#ifndef _MHDPARAMS_H_
#define _MHDPARAMS_H_

class IniDocument;
class MHDParams
{   
    public:
        double gamma;
        void init(IniDocument& ini_doc);
        void write_mhdhelp(const char* filename);
};
extern MHDParams mhdParams;

#endif
