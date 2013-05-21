#ifndef RayIntersectableEllipse_h
#define RayIntersectableEllipse_h

#include <cmath>

#include <DGtal/kernel/NumberTraits.h>
#include <DGtal/base/BasicFunctors.h>


/**
 * Class implementing a circle that is 'ray intersectable', 
 * ie. the intersection between the circle and 
 * a ray emanating from a given point along a given direction
 * is computable.
 * 
 * @tparam Point a model of point/vector.  
 * @tparam Parameter a model of integer for the circle parameters
 *
 * This class is a model of point functor and of ray intersectable shape
 */
template <typename TPoint, typename TParameter = long long int>
class RayIntersectableEllipse 
{
  public: 
    /////////////////////// inner types /////////////////
    typedef TPoint Point; 
    typedef TPoint Vector; 
    typedef typename TPoint::Coordinate Coordinate; 
    typedef TParameter Integer; 
    typedef Integer Value; 

  private: 
    /////////////////////// members /////////////////////
    /**
     * Circle parameters. The circle is defined
     * as the set {(x,y) | ax + by + c(x^2 + y^2) + d = 0}
     * The center is equal to ( -a/(2c) , -b/(2c) ) and the 
     * radius is equal to sqrt( (a^2 + b^2 - 4cd) / 4c^2 )
     */
    Integer myA, myB, myC, myD, myE, myF;  

  public:
    ///////////////////// standard services /////////////
    /**
     * Default constructor
     */
    RayIntersectableEllipse(): myA(0), myB(0), myC(0), myD(0), myE(0), myF(0) {}

    /**
     * Standard constructor
     * @param a a-parameter
     * @param b b-parameter
     * @param c c-parameter
     * @param d d-parameter
     */
    RayIntersectableEllipse(const Integer& a, const Integer& b, const Integer& c, 
      const Integer& d, const Integer& e, const Integer& f)
      : myA(a), myB(b), myC(c), myD(d), myE(e), myF(f) {}


    /**
     * Copy constructor
     * @param p other object to copy
     */
    RayIntersectableEllipse(const RayIntersectableEllipse& p): 
      myA(p.a()), myB(p.b()), myC(p.c()), myD(p.d()), myE(p.e()), myF(p.f()) {}

    /**
     * Assignement operator
     * @param p other object to copy
     * @return reference on *this
     */
    RayIntersectableEllipse& operator=(const RayIntersectableEllipse& p) 
    {
      if (p != *this)
      {
        myA = p.a(); 
        myB = p.b(); 
        myC = p.c(); 
        myD = p.d();
        myE = p.e();
        myF = p.f();        
      }
      return *this; 
    }

    /**
     * Default destructor
     */
    ~RayIntersectableEllipse() {}

    ///////////////////// read access ///////////////////
    /**
     * a-parameter accessor
     * @return a-parameter. 
     */
    Integer a() const { return myA; }

    /**
     * b-parameter accessor
     * @return b-parameter. 
     */
    Integer b() const { return myB; }

    /**
     * c-parameter accessor
     * @return c-parameter. 
     */
    Integer c() const { return myC; }

    /**
     * d-parameter accessor
     * @return d-parameter. 
     */
    Integer d() const { return myD; }

    /**
     * e-parameter accessor
     * @return e-parameter. 
     */
    Integer e() const { return myE; }

    /**
     * f-parameter accessor
     * @return f-parameter. 
     */
    Integer f() const { return myF; }




  public: 
    ///////////////////// main methods ///////////////////
    /**
     * Function operator
     * @param aPoint any point
     * @return 0 if @a aPoint is on the ellipse,
     * a value <0 if @a aPoint is outside
     * a value >0 if @a aPoint is inside
     */
    Value operator()(const Point& aPoint) const
    {
      Integer x = aPoint[0]; 
      Integer y = aPoint[1];
      Integer xy = x*y;
      Integer xx = x*x;
      Integer yy = y*y;

      return (myA*xx + myB*xy + myC*yy + myD*x + myE*y + myF); 
    }

  private:
    /**
     * Retrieves integer @a q such that 
     * the point @a aPoint + @a q * @a aDir 
     * and the points @a aPoint + (@a q+1) * @a aDir
     * lie in either side of the circle. 
     * @param aPoint source of the ray casting
     * @param aDir direction of the ray casting
     * @param qmax maximal integer bounding the 
     * dichotomic search
     * @return integer @a q
     */
    template <typename Predicate>
      Coordinate dichotomicSearch(const Point& aPoint, const Vector& aDir, 
          const Coordinate& qmax, const Predicate& aPredicate) const
      {    

        ASSERT( aPredicate( (*this)(aPoint)) ); 
        ASSERT( (!aPredicate( (*this)(aPoint + qmax * aDir))) ); 

        // init search bounds
        Coordinate qStart = 0;
        Coordinate qStop  = qmax;
        Point pStart = aPoint; 
        Point pStop = aPoint + qmax * aDir; 

        Coordinate qMid;
        Point pMid, pMid1; 

        // while not yet located 
        bool isNotFound = true; 
        while( (qStop > qStart) && (isNotFound) ) 
        {
          // middle between qStart and qStop
          qMid = (qStart + qStop)/2;
          pMid = aPoint + qMid * aDir;
          pMid1 = aPoint + (qMid+1) * aDir;

          // inclusion tests
          if ( aPredicate( (*this)(pMid)) )
          { 
            if (!aPredicate( (*this)(pMid1)) )
            { //found!
              isNotFound = false; 
              qStart = qMid; 
            }
            else //search in the upper range
              qStart = qMid + 1;
          }
          else //search in the lower range
            qStop = qMid;
        }
        // return the integer 
        return(qStart);
      }

  public: 
    /**
     * Intersection between the circle and the digital ray
     *
     * @param aStartingPoint point from which the ray is emanating
     * @param aDirection ray direction
     * @param returned maximal integer such that
     * @a aStartingPoint + @a aQuotient * @a aDirection lies close
     * (but on the same side) to the circle
     * @param returned point equal to @a aStartingPoint + @a aQuotient * @a aDirection
     *
     * @return 'true' if the ray and the circle instercest, 'false' otherwise 
     */
    bool dray(const Point& aStartingPoint, const Vector& aDirection, 
        Coordinate& aQuotient, Point& aClosest) const 
    {

//We look for the intersection between  aS + aQ * aD and the ellipse
      Integer aD0 = aDirection[0];
      Integer aD1 = aDirection[1];

      Integer aS0 = aStartingPoint[0];
      Integer aS1 = aStartingPoint[1];    


      Integer aEq = myA*aD0*aD0 + myB*aD0*aD1 + myC*aD1*aD1;
      Integer bEq = 2*myA*aS0*aD0 + myB*(aS0*aD1 + aS1*aD0) + 2*myC*aS1*aD1 + myD*aD0 + myE*aD1;
      Integer cEq = myA*aS0*aS0 + myB*aS0*aS1 + myC*aS1*aS1 + myD*aS0 + myE*aS1 + myF;

      Integer Delta = bEq*bEq - 4*aEq*cEq;

      //Delta = 0 : 1 solution
      // aS+aQ*aD is a tangent to the circle
      if ( (Delta == 0) && (-bEq/(2*aEq) >= 0) )
      {
        aQuotient = DGtal::NumberTraits<Integer>::castToInt64_t(-bEq /(2*aEq));
        aClosest = aStartingPoint + aQuotient*aDirection;
        return true;
      }
      // Delta > 0 : 2 solutions
      else if (Delta > 0)
      {
        if ( cEq == 0)
        { // aS lies on the circle
          aQuotient = 0;
          aClosest = aStartingPoint;
          return true;
        }
        else 
        {  
          if ( cEq < 0)
          { // aS is strictly outside the circle
            if ( (-bEq/aEq) < 0  ){return false;}
            else 
            { // If positive solutions, 
              // We pick the smallest one by dichotomic search
              //this predicate is true iff the input value is <= 0
              using namespace DGtal; 
              Thresholder<Integer, true, true> predicate(0);

              Coordinate max = DGtal::NumberTraits<Integer>::castToInt64_t(-bEq/(2*aEq)); 


              //inclusion tests around the minimum (-bEq/(2*aEq))
              Point pMax = aStartingPoint + max * aDirection; 
              Integer pMaxValue = (*this)(pMax);  
              Point pMax1 = aStartingPoint + (max+1) * aDirection;
              Integer pMaxValue1 = (*this)(pMax1);  
              if ( pMaxValue < 0 ) 
              { //if pMax is strictly outside the shape
                if ( pMaxValue1 < 0 )
                { //if pMax1 is outside the shape
                  return false; //there is no integral solution
                }
                else 
                {
                  if ( pMaxValue1 == 0 ) //if on the boundary
                    aQuotient = max + 1; 
                  else        //if strictly inside the circle
                    aQuotient = dichotomicSearch(aStartingPoint, aDirection, (max+1), predicate);
                }
              }
              else
              {
                if ( pMaxValue == 0 ) //if on the boundary
                  aQuotient = max; 
                else        //if strictly inside the circle
                  aQuotient = dichotomicSearch(aStartingPoint, aDirection, max, predicate);
              }

              aClosest = aStartingPoint + aQuotient*aDirection;
              return true;
            }
          }
          else
          { // aS is strictly inside the circle
            // We pick the positive solution by dichotomic search
            //this predicate is true iff the input value is >= 0
            /*using namespace DGtal;
            Thresholder<Integer, false, true> predicate(0); 
            //a trivial upper bound is the circle diameter
            Coordinate diameter = 2 * (Coordinate) std::ceil( getRadius() ) + 1; 
            aQuotient = dichotomicSearch(aStartingPoint, aDirection, diameter, predicate); 
            aClosest = aStartingPoint + aQuotient*aDirection;*/
            
            aQuotient = 0;
            while ((*this)(aStartingPoint + aQuotient*aDirection) > 0)
            {
              aQuotient++;
            }
            aQuotient--; 
            aClosest = aStartingPoint + aQuotient*aDirection;
            
            return true;
          }
        }
      }
      // Delta < 0 and cases where there is no solution
      else {return false;}
    }

    /**
     * Returns a vertex of the convex hull of the digital points
     * lying inside the circle, such that it has a minimal y-coordinate, 
     * and among the digital points of minimal y-coordinate, it has
     * a maximal x-coordinate. 
     * @return vertex of the convex hull. 
     */
    Point getConvexHullVertex() const 
    { 
      return(Point(0,0));
    }

}; 
#endif

