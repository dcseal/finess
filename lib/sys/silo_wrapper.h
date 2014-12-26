#ifndef FINESS_SILO_H
#define FINESS_SILO_H

#include <silo.h>
#include <stdexcept>
#include <string>


namespace finess {
    namespace silo {
        using std::runtime_error;
        using std::string;
        class Exception: runtime_error {
            public:
                explicit Exception(const string& what_arg) :
                    runtime_error(what_arg)   
            { }
        };
        class SetCompression {
            private:
                string old_options;
            public:
                SetCompression(string options){
                    const char* s = DBGetCompression();
                    old_options = (s == NULL) ? "" : s;
                    DBSetCompression(options.c_str());
                }
                ~SetCompression(){
                    DBSetCompression(old_options == "" ? NULL : old_options.c_str());
                }
        };
        class File {
            private:
                DBfile *dbfile;
            public:
                File(string filename, string comment = "")
                {
                    dbfile = NULL;
                    dbfile = DBCreate(filename.c_str(), DB_CLOBBER, DB_LOCAL, comment.c_str(), DB_HDF5);
                    if(dbfile == NULL)
                        throw Exception( string("Could not create Silo file: ") + filename );
                }
                ~File(){
                    DBClose(dbfile);
                }
                DBfile *get_ptr(){
                    return dbfile;
                }
        };
        class TimeOptions {
            private:
                DBoptlist *optlist;
            public:
                TimeOptions(int nframe, double time){
                    optlist = DBMakeOptlist(2);
                    if(optlist == NULL)
                        throw Exception(string("Could not create Silo optlist."));
                    DBAddOption(optlist, DBOPT_CYCLE, &nframe);
                    DBAddOption(optlist, DBOPT_DTIME, &time);
                }
                ~TimeOptions(){
                    DBFreeOptlist(optlist);
                }
                DBoptlist* get_ptr() const{
                    return optlist;
                }
        };
    }
}

#endif
