#ifndef tensors1d_h
#define tensors1d_h

#ifndef NDIMS
#include <dimdefs.h> // for NDIMS (number of spatial dimensions)
#endif
#if (!defined(NDIMS) || !(NDIMS>=0))
#error "NDIMS must be defined"
#endif

// What's the difference between DEEP and DIMS copy? (-DS)
namespace CopyMode
{
    enum Enum
    {
        DEEP = 1,
        DIMS = 2,
    };
}



/**
@brief Wrapper class for a flat array of integer type.

This class is to be used only as a base class to other classes.

- The <tt>iTensor</tt> classes are used only in #fprint_tensor(...), #str_into_tensor(...);
- Uses a heavily compiler-dependent mechanism to do bounds-checking.  @sa CHECK_BOUNDS.
*/
class iTensorBase
{
    protected:
        void init();
        ///@brief Default constructor.  Protected so that only derived class can call.
        ///
        ///Actually, not called anywhere.
        iTensorBase(){};
        iTensorBase(const iTensorBase& in);
    public:
        iTensorBase(int size_in):size(size_in){init();}
        ~iTensorBase();
        void setall(int);
        const int numel() const { return size; }
#ifdef CHECK_BOUNDS
        const int& vget(int k) const;
        int& vfetch(int k);
        void vset(int k, int value);
#else
        const int& vget(int k) const
        {
            return vec[k];
        }
        int& vfetch(int k) {
            return vec[k];
        }
        void vset(int k, int value){
            vec[k]=value;
        }
#endif
    protected:
        int* vec;
        int size;
};


///@brief A 1d integer array with designated starting index and size.
///
///- Never used anywhere.
class iTensor1d : public iTensorBase
{
    // data
    protected:
        int b1;
        // methods
    protected:
        iTensor1d(){};
    public:
        int getidx(int n1) const
        {
            int k = n1-b1;
            return k;
        }
        int getsize() const { return size; }
    public:
        // constructor takes size and initial index in each dimension
        iTensor1d( int s1i, int b1i ) : b1(b1i) { size=s1i; init(); }
        iTensor1d(const iTensor1d& in) : iTensorBase(in), b1(in.b1) {}

        const int& get(int n1) const
        { return vget(n1-b1); }
        void set(int n1, int value)
        { return vset(n1-b1, value); }
};

///@brief A 1d integer array, with starting index=1.
///
///- Never used anywhere.
class iTensor1 : public iTensor1d
{
    public:
        iTensor1(int s1i) :
            iTensor1d(s1i,1) { }
};

// === major section: dtensors (arrays of double) ===

///@brief A flat array of <tt>double</tt> type.  Starting index=0.
///
///This class is to be used only as a base class to other tensor classes.
///
///- Bounds-check uses a heavily compiler-dependent mechanism.  See #CHECK_BOUNDS.
///  As a consequence, #vget and #vfetch in the current document do not necessarily correspond to
///  where they are actually defined (after the preprocessing, and compiler-dependent behavivor of linking).
///@todo Clean up bounds-check.
class dTensorBase
{
    private: // disabled
        ///@todo Replace with <tt>=delete</tt> from C++11 Standard.
        dTensorBase& operator=(const dTensorBase& in);
    protected:
        void init();
        ///@brief Protected default constructor.
        ///
        ///- Never called anywhere.
        dTensorBase(){};
        dTensorBase(const dTensorBase& in, CopyMode::Enum copyMode=CopyMode::DEEP);
        void copyfrom(const dTensorBase& in);
    public:
        ///@brief Set #size to <tt>size_in</tt> and calls #init(...).
        dTensorBase(int size_in):size(size_in){init();}
        ~dTensorBase();
        bool check();
        void setall(double);
        ///@brief Returns number of elements.
        const int numel() const { return size; }
#ifdef CHECK_BOUNDS
        const double& vget(int k) const;
        double& vfetch(int k);
        void vset(int k, double value);
#else
        ///@note Two definitions
        ///@todo Clean up re-definitions.
        const double& vget(int k) const { return vec[k]; }
        ///@note Two definitions
        ///@todo Clean up re-definitions.
        double& vfetch(int k) {
            return vec[k];
#ifdef CHECK_INIT
            if(!(vec[k]==vec[k])) eprintf("vec[%d]=%24.16e",k,vec[k]);
#endif
        }
        void vset(int k, double value){vec[k]=value;}
#endif
    protected:
        ///@brief Pointer holding starting address of array.
        double* vec;
        ///@brief Number of elements in the array.
        int size;
};


///@brief One-dimensional array with arbitrary beginning index and size.
///
///- The indexing of the array by {desired beginning index, ..., desired beginning index + number of elements - 1}
///  will be referred to as "high-level" indexing.
///  The indexing by {0, ..., number of elements - 1} will be referred to as "low-level" indexing.
///- #copyfrom(...) from #dTensorBase (which was protected) is exposed to public, with additional guard (assert equal) on beginning index.
class dTensor1d : public dTensorBase
{
    // data
    protected:
        ///@brief Beginning index.
        int b1;
        // methods
    private: // disabled
        ///@brief Disabled copy assignment operator.
        ///@todo Replace with <tt> =delete </tt> in C++11 Standard.
        dTensor1d& operator=(const dTensor1d& in);
    protected:
        ///@brief Protected default construtor.
        dTensor1d(){};
    public:
        ///@brief Returns low-level index corresponding to high-level index <tt>n1</tt>.
        int getidx(int n1) const
        {
            int k = n1-b1;
            return k;
        }
        ///@brief Returns the number of elements.
        int getsize() const { return size; }
    public:
        // constructor takes size and initial index in each dimension

        ///@brief Constructs an array of <tt>s1i</tt> elements, with beginning index <tt>b1i</tt>.
        dTensor1d( int s1i, int b1i ) : b1(b1i) { size=s1i; init(); }
        ///@brief Copy constructor.
        ///
        ///- Copies both the underlying flat array (by calling copy constructor of #dTensorBase)
        ///  and the beginning index.
        dTensor1d(const dTensor1d& in) : dTensorBase(in), b1(in.b1) {}
        void copyfrom(const dTensor1d& in);

        ///@brief Returns (<tt>const</tt> reference to) <tt>n1</tt>-th element in the high-level indexing.
        const double& get(int n1) const
        { return vget(n1-b1); }
        ///@brief Returns reference to <tt>n1</tt>-th element in the high-level indexing.
        double& fetch(int n1)
        { return vfetch(n1-b1); }
        ///@brief Sets <tt>n1</tt>-th element (in the high-level indexing) to <tt>value</tt>.
        void set(int n1, double value)
        { return vset(n1-b1, value); }
};



///@brief One-dimensional array with beginning index=1.
///
///Difference from #dTensor1d:
///- Constructor takes the number of elements as the only argument;
///- #copyfrom(...) is re-defined to accept a source of type #dTensor1.
class dTensor1 : public dTensor1d
{
    public:
        ///@brief Constructs a one-dimensional array with beginning index=1, and number of elements=<tt>s1i</tt>
        dTensor1(int s1i) :
            dTensor1d(s1i,1) { }
        ///@brief Copies from <tt>in</tt>.
        ///
        ///- Implemented by invoking #dTensor1d::copyfrom(...).  Only difference: argument type is more specialized.
        void copyfrom(const dTensor1& in){ dTensor1d::copyfrom(in); }
    private: // disabled
        ///@brief Copy assignment operator that is disabled (to the public).
        ///
        ///@todo  Clarify:  is it used inside the class?  If not, why not an empty body?
        dTensor1& operator=(const dTensor1& in){ copyfrom(in); return *this; }
};



///@brief One-dimensional array with boundary cells.
///
///Difference from #dTensor1d:
///- Arguments to constructor is interpreted differently;
///- Additional internal variables: #S1, #mbc, and #ndims;
///- Additional asserts in #copyfrom(...);
///- Added #getmbc(...) method.
class dTensorBC1 : public dTensor1d
{
    private: // disabled
        ///@brief Disabled copy assignment operator.
        ///@todo Replace with <tt>=delete</tt>.
        dTensorBC1& operator=(const dTensorBC1& in);
    public:
        ///@brief Constructs 
        dTensorBC1(int s1i,
                int mbc, int ndims=1);
        void copyfrom(const dTensorBC1& in);
        int getmbc() const {return mbc;}
    private:
        ///@brief Number of cells that are not boundary cells.
        int S1;
        ///@brief Number of boundary cells that are padded on the two ends of each dimension that needs such padding.
        int mbc;
        ///@brief The number of dimensions that need paddings on two ends for boundary cells.
        int ndims;
};

// methods for tensors
//The following three functions are defined in lib/dog_str.cpp
int count_fields(const char* str_in, char field_sep);
bool str_into_tensor(const char* str, iTensorBase& t, char field_sep);
bool str_into_tensor(const char* str, dTensorBase& t, char field_sep);
#endif
