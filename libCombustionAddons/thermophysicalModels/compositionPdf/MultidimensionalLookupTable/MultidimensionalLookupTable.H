#ifndef MULTIDIMENSIONALLOOKUPTABLE_H
#define MULTIDIMENSIONALLOOKUPTABLE_H

#include "Istream.H"
#include "Ostream.H"
#include "Field.H"

#ifndef SMALL
#define SMALL 1e-10
#endif


double Round(double Zahl, int Stellen=0);
     

namespace Foam
{

    template<class T>
    class MultidimensionalLookupTable;



    template<class T>
    Ostream& operator<<
        (
            Ostream& f, 
            const MultidimensionalLookupTable<T>& t
        );



    template<class T>
    class MultidimensionalLookupTable
    {

    protected:
        label dims;
        scalar *start;
        scalar *delta;
        label *num;
        T* elements;

    public:


      static autoPtr<MultidimensionalLookupTable<T> > New(Istream& is)
      {
        return autoPtr<MultidimensionalLookupTable<T> >
            (new MultidimensionalLookupTable<T>(is));
      }


      virtual autoPtr<MultidimensionalLookupTable<T> > clone() const
      {
        return autoPtr<MultidimensionalLookupTable<T> >
            (new MultidimensionalLookupTable<T>(*this));
      }


      inline label dimensionality() const { return dims; };

      MultidimensionalLookupTable
          (
           const MultidimensionalLookupTable& o
          )
        : dims(o.dims)
      {
        start=new scalar[dims];
        delta=new scalar[dims];
        num=new label[dims];

        memcpy(start, o.start, sizeof(scalar)*dims);
        memcpy(delta, o.delta, sizeof(scalar)*dims);
        memcpy(num, o.num, sizeof(scalar)*dims);
        
        int size=num[0];
        for (int i=1;i<dims;i++) size*=num[i];
        elements=new T[size];
        for (int i=1;i<size;i++)
          elements[i]=o.elements[i];
      }

        MultidimensionalLookupTable
        (
            const label dim,
            const label *n, 
            const scalar* s, 
            const scalar* d
        )
            : dims(dim)
            {

                start=new scalar[dims];
                delta=new scalar[dims];
                num=new label[dims];

                // read header
                for (int i=0;i<dims;i++) num[i]=n[i];
                for (int i=0;i<dims;i++) start[i]=s[i];
                for (int i=0;i<dims;i++) delta[i]=d[i];
                
                // alloc
                int size=num[0];
                for (int i=1;i<dims;i++) size*=num[i];
                elements=new T[size];
            }
  
        MultidimensionalLookupTable
        (
            Istream& s
        )
            {

                token beginToken(s); // leading "("

                s>>dims;
                start=new scalar[dims];
                delta=new scalar[dims];
                num=new label[dims];
    
                // read header
                for (int i=0;i<dims;i++) s>>num[i];
                for (int i=0;i<dims;i++) s>>start[i];
                for (int i=0;i<dims;i++) s>>delta[i];

                // alloc
                int size=num[0];
                for (int i=1;i<dims;i++) size*=num[i];
                elements=new T[size];

                // read data
                T elem;
                int idx[dims];
                memset(idx, 0, sizeof(idx));
                for (label i=0;i<size;i++)
                {
                    s>>elem;
                    if (s.fail()) break;
                    operator[](idx)=elem;
                    //Info<<"saved in: "<<idx[0]<<" "<<idx[1]<<endl;
                    incIdx(idx);
                }
                
                token endToken(s); // trailing ")"

            }

        virtual ~MultidimensionalLookupTable()
            {
                delete start;
                delete delta;
                delete num;
	        //delete elements;
            }

        inline label nElements(label i) const { return num[i]; }
        inline scalar startValue(label i) const { return start[i]; }
        inline scalar deltaValue(label i) const { return delta[i]; }
        inline scalar maxValue(label i) const 
            { return start[i] + scalar(num[i]-1)*delta[i]; }

        inline void rescale(label i, scalar fac)
            { delta[i]*=fac; }

        inline scalar valueAt(label dim, label i) const
            { return start[dim]+scalar(i)*delta[dim]; }

        inline label indexAt(label dim, scalar value) const
            { scalar i=(value-start[dim])/delta[dim]; Info<<i<<endl; return label(Round(i));}

        inline T& access(const int *idx)
            {
                int j=0, m=1;
                for (int i=0;i<dims;i++)
                {
                    //int bi=max(min(idx[i],num[i]-1),0);
                    int bi=idx[i];
                    if (bi>num[i]-1)
                    {
                     FatalErrorIn("MultidimensionalLookupTable::access()")
                      << "Index of dimension " << i
                      << " ("<<idx[i]<<") exceeds upper bound " <<num[i]-1 <<endl
                      << abort(FatalError);
                    }
                    if (bi<0)
                    {
                     FatalErrorIn("MultidimensionalLookupTable::access()")
                      << "Index of dimension " << i
                      << " ("<<idx[i]<<") exceeds lower bound 0" <<endl
                      << abort(FatalError);
                    }
                    j+=bi*m;
                    m*=num[i];
                }
                return elements[j];
            }
        
        inline const T& const_access(const int *idx) const
            {
                int j=0, m=1;
                for (int i=0;i<dims;i++)
                {
                    int bi=max(min(idx[i],num[i]-1),0);
                    j+=bi*m;
                    m*=num[i];
                }
                return elements[j];
            }

        inline T& operator[](const int* idx)
            {
                return access(idx);
            }

        inline const T& operator[](const int *idx) const
            { 
                return const_access(idx); 
            }
        
        bool incIdx(int *idx) const
            {
                for (int i=0;i<dims;i++)
                {
                    if (idx[i]<num[i]-1)
                    {
                        idx[i]++;
                        return true;
                    }
                    else
                    {
                        idx[i]=0;
                    }
                }
                return false;
            }
  
  
        void write(Ostream& f) const
            {
                
                f<<token::BEGIN_LIST<<endl;

                f<<dims<<token::SPACE;

                // write header
                for (int i=0;i<dims;i++) f<<num[i]<<token::SPACE;
                f<<endl;
                for (int i=0;i<dims;i++) f<<start[i]<<token::SPACE;
                f<<endl;
                for (int i=0;i<dims;i++) f<<delta[i]<<token::SPACE;
                f<<endl;
                
                int size=num[0];
                for (int i=1;i<dims;i++) size*=num[i];
                
                // write data
                int idx[dims];
                memset(idx, 0, sizeof(idx));
                for (label i=0;i<size;i++)
                {
                    f<<const_access(idx)<<endl;
                    incIdx(idx);
                }
                
                f<<token::END_LIST<<endl;  
                
            }

        List<scalar>  compositionTableGradient
        (
            const List<scalar>& pp, 
            const label direction=0
        ) const
            {
                List<scalar> p(pp);
                
                for (label i=0;i<dims;i++)
                {
                    if (p[i]<start[i])
                    {
                        p[i]=start[i];
                    }
                    if (p[i]>num[i]*delta[i])
                    {
                        p[i]=num[i]*delta[i];
                    }
                }


                int idx0[dims];
                
                for (int i=0;i<dims;i++)
                {
                    double pos=p[i];
                    idx0[i]=int(floor((pos-start[i])/delta[i]));
                }
                
                Field<scalar> ret=operator[](idx0).composition();
                if ( idx0[direction] < num[direction] )
                {
                    idx0[direction]++;
                    ret*=-1.0;
                    ret+=operator[](idx0).composition();
                } 
                else
                {
                    idx0[direction]--;
                    ret-=operator[](idx0).composition();
                }
                ret/=delta[direction];
                
                return ret;
            }
 
#define	WARNINGTHRESHOLD 1e-3
 
        T lookup(const List<scalar>& pp) const
            {
                List<scalar> p(pp);
                
                for (label i=0;i<dims;i++)
                {
                    if ( start[i] > p[i])
                    {
/*
                        if (start[i]-p[i]>WARNINGTHRESHOLD)
                        Info
                            //    FatalErrorIn("MultidimensionalLookupTable::lookup()")
                            <<"Table Lookup: value ("
                                <<i<<") out of lower bound ("<<start[i]<<"): "<<p[i]<<endl;
                            //  <<abort(FatalError);
*/
                        p[i]=start[i];
                    }
                    if ( p[i] > start[i]+(num[i]-1)*delta[i] )
                    {
/*
                        if (p[i]-(start[i]+(num[i]-1)*delta[i])>WARNINGTHRESHOLD)
                        Info<<"Table Lookup: value ("
                            <<i<<") out of upper bound ("<<start[i]+(num[i]-1)*delta[i]<<"): "<<p[i]<<endl;
*/
                        p[i]=start[i]+(num[i]-1)*delta[i];
                    }
                }

                int idx0[dims];
                double d[2][dims];
                
                for (int i=0;i<dims;i++)
                {
                    double pos=p[i];
                    idx0[i]=int(floor((pos-start[i])/delta[i]));
                    d[0][i]=pos-(start[i]+idx0[i]*delta[i]);
                    d[1][i]=start[i]+(idx0[i]+1)*delta[i]-pos;
                }
                
                T ret;
                
                    #define OT(x) (x==0?1:0)
                
                double sumfac=0.0;
		bool assigned=false;
                for (int c=0;c<(1<<dims);c++)
                {
                    int idxl[dims], idx[dims];
                    for (int b=0;b<dims;b++)
                    {
                        idxl[b]=((c&(1<<b))>0);
                        idx[b]=idx0[b]+idxl[b];
                    }

                    double fac=1.0;
                    for (int b=0;b<dims;b++)
                    {
                        fac*=d[OT(idxl[b])][b];
                    }
                    sumfac+=fac;
		    if (assigned) 
		      ret+=fac*(operator[](idx));
		    else 
		      { ret=fac*(operator[](idx)); assigned=true; }
                }
                ret/=sumfac;
                return ret;
            }


        friend
        Ostream& operator<< <T> 
                            (
                                Ostream& f, 
                                const MultidimensionalLookupTable<T>& t
                            );
        
    };


    template<class T>
    Ostream& operator<<
        (
            Ostream& f, 
            const MultidimensionalLookupTable<T>& t
        )
    {        
        t.write(f);        
        return f;
    }
      

}; // Namespace Foam

#endif
