/**
 * 
 */
package de.hzdr.jgm.cgeo.gmap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import de.hzdr.jgm.cgeo.Base;

//import cern.colt.matrix.impl.DenseDoubleMatrix1D;
//import cern.colt.matrix.impl.DenseDoubleMatrix2D;
//import cern.colt.matrix.linalg.Algebra;

/**
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class Sphere
{

   private double[] center =
   { 0, 0, 0 };
   private double   radius = 0;

   public Sphere()
   {
      super();
   }

   public Sphere(double[] center, double radius)
   {
      super();
      this.center = center;
      this.radius = radius;
   }

   public Sphere(double[] point)
   {
      super();
      center = point;
      radius = 0.0;
   }

   public Sphere(double[] point1, double[] point2)
   {
      super();
      center = new double[point1.length];
      for (int i = 0; i < point1.length; i++)
      {
         center[i] = (point1[i] + point2[i]) * 0.5;
         radius = radius + point1[i] * point2[i];
      }
      radius = Math.sqrt(radius);

   }

   /**
    * @return the center
    */
   public double[] getCenter()
   {
      return center;
   }

   /**
    * @param center
    *           the center to set
    */
   public void setCenter(double[] center)
   {
      this.center = center;
   }

   /**
    * @return the radius
    */
   public double getRadius()
   {
      return radius;
   }

   /**
    * @param radius
    *           the radius to set
    */
   public void setRadius(double radius)
   {
      this.radius = radius;
   }

   /**
    * point in sphere test
    * 
    * @param point
    * @return
    */
   public boolean isIn(double[] point)
   {
      if (point.length == center.length)
      {
         return (Base.distanceSq(center, point) - radius * radius) < GMap.TOLERANCE_EPS;
      }
      else
         throw new IndexOutOfBoundsException("Inconsistent dimensionality!");
   }

   /**
    * point in sphere test
    * 
    * @param point
    * @param tolerance
    * @return
    */
   public boolean isIn(double[] point, double tolerance)
   {
      if (point.length == center.length)
      {
         return (Base.distanceSq(center, point) - radius * radius) < tolerance;
         // return (Base.distanceSq(center, point))<((radius*radius)+tolerance);
      }
      else
         throw new IndexOutOfBoundsException("Inconsistent dimensionality!");
   }

   /**
    * Creates a minimum-circumsphere around a set of points
    * 
    * @param points
    *           set of points
    * @return returns the list of points located ON the sphere
    */
   public List<double[]> defineBy(List<double[]> points)
   {
      boolean found_sph = false;

      int m_dim = points.get(0).length;
      List<double[]> ret = new ArrayList<double[]>();

      Sphere sph = null;
      if (points.size() == 2)
      {
         sph = new Sphere(points.get(0), points.get(1));
         ret.add(points.get(0));
         ret.add(points.get(1));
         found_sph = true;
      }
      else
      {
         // level 2 - find max. distance between 2 points
         int level = 2;
         ArrayList<List<Integer>> combinations = Base.getAllIndexCombinations(points.size(), level);
         double maxDist = 0;
         int id = 0;

         for (int i = 0; i < combinations.size(); i += 1)
         {
            double d = Base.distanceSq(points.get(
                  combinations.get(i).get(0)),
                  points.get(combinations.get(i).get(1)));

            if (d > maxDist)
            {
               maxDist = d;
               id = i;
            }
         }

         sph = new Sphere(points.get(combinations.get(id).get(0)),
               points.get(combinations.get(id).get(1)));

         found_sph = sph.containsCoordinateList(points);

         level = 3;
         while (!found_sph && level <= m_dim + 1)
         {
            combinations = Base.getAllIndexCombinations(points.size(), level);

            double minRad = -1;
            for (int c = 0; c < combinations.size(); c += 1)
            {
               ArrayList<double[]> pts = new ArrayList<double[]>();
               for (int i = 0; i < combinations.get(c).size(); i += 1)
               {
                  pts.add(points.get(combinations.get(c).get(i)));
               }

               Sphere sph2 = Sphere.getCircumscribedSphere(pts, true);

               if (sph2 != null && (minRad == -1 || sph2.getRadius() < minRad))
               {
                  if (sph.containsCoordinateList(points))
                  {
                     ret = points;
                     sph = sph2;
                     found_sph = true;
                     minRad = sph2.getRadius();
                  }
               }
            }
            level++;
         }
      }
      center = sph.center;
      radius = sph.radius;
      return ret;
   }

   public boolean containsCoordinateList(List<double[]> points)
   {
      for (int p = 0; p < points.size(); p += 1)
      {
         if (!this.isIn(points.get(p)))
            return false;
      }
      return true;
   }

   /**
    * Creates a minimum-circumsphere around a set of points and based on a single point located on its boundary
    * 
    * @param points
    *           set of points
    * @param p
    *           point on sphere boundary
    * 
    * @return returns the list of points located ON the sphere
    */
   public List<double[]> defineBy(List<double[]> points, double[] p)
   {
      return null;
   }

   public String print()
   {
      String ret = this + ":\n\tcenter: " + Arrays.toString(center) + ":\n\tradius: " + radius;
      return ret;
   }

   /**
    * Creats a N-Sphere as circumsphere of a set of n-points using Colt.linalg - functionality - all points are located
    * on the sphere
    * 
    * @param points
    * @return new n-Sphere defined by the given points or null, when no sphere can be obtained
    */
   public static Sphere getCircumscribedSphere(List<double[]> points, boolean useColt)
   {
      int k = points.size() - 1;

      if (points.size() > 0)
      {
         int m_dim = points.get(0).length;
         double rds = 0;
         double[] cntr = new double[m_dim];

         if (points.size() == 1)
         {
            // trivial case I: only one point in the center, radius = 0 (for all dimensions)
            cntr = points.get(0);
         }
         else
            if (points.size() == 2)
            {
               double[] point1 = points.get(0);
               double[] point2 = points.get(1);
               for (int i = 0; i < m_dim; i++)
               {
                  cntr[i] = (point1[i] + point2[i]) * 0.5;
                  rds = rds + point1[i] * point2[i];
               }
               rds = Math.sqrt(rds);
            }
            else
               if (k == m_dim)
               {
                  // standard case: formulas from Cheng, Hu & Martin (2006): ON THE SMALLEST ENCLOSING BALLS
                  double[][] A = new double[m_dim][m_dim]; // matrix A, should be invertible, or the given
                  // // points are not located all on the sphere
                  double[] B = new double[m_dim];

                  for (int i = 0; i < k; i += 1)
                  {
                     double bj = 0;
                     for (int j = 0; j < m_dim; j += 1)
                     {
                        A[i][j] = points.get(i + 1)[j] - points.get(i)[j];
                        bj = bj + 0.5 * (points.get(i + 1)[j] * points.get(i + 1)[j]
                              - points.get(i)[j] * points.get(i)[j]);
                     }
                     B[i] = bj;
                  }

//                  if (useColt)
//                  {
//                     DenseDoubleMatrix2D Addm = new DenseDoubleMatrix2D(A);
//                     Algebra a = new Algebra();
//                     if (a.det(Addm) != 0)
//                     {
//                        DenseDoubleMatrix1D Bddm = new DenseDoubleMatrix1D(B);
//                        cntr = a.mult(a.inverse(Addm), Bddm).toArray();
//                        rds = Base.distance(cntr, points.get(0));
//                     }
//                     else
//                     {
//                        return null;
//                     }
//                  }
//                  else
                  {
                     if (m_dim == 2)
                        cntr = Base.solveExplicit2x2(A, B);
                     else
                        if (m_dim == 3)
                           cntr = Base.solveExplicit3x3(A, B);
                        else
                           cntr = Base.solve_Cramer(A, B, GMap.TOLERANCE_EPS);

                     if (cntr != null)
                     {
                        rds = Base.distance(cntr, points.get(0));
                     }
                     else
                        return null;
                  }
               }
               else
                  if (k < m_dim)
                  {
                     // non trivial cases form Cheng, Hu & Martin (2006): ON THE SMALLEST ENCLOSING BALLS
                     double[][] A = new double[m_dim][m_dim]; // matrix A, should be invertible, or the
                     // given
                     // // points are not located all on the sphere
                     double[] B = new double[m_dim];
                     // first part, same as before, but with to less rows
                     for (int i = 0; i < k; i += 1)
                     {
                        double bj = 0;
                        for (int j = 0; j < m_dim; j += 1)
                        {
                           A[i][j] = points.get(i + 1)[j] - points.get(i)[j];
                           bj = bj + 0.5 * (points.get(i + 1)[j] * points.get(i + 1)[j]
                                 - points.get(i)[j] * points.get(i)[j]);
                        }
                        B[i] = bj;
                     }

                     // second part
//                     if (useColt)
//                     {
//                        Algebra a = new Algebra();
//                        DenseDoubleMatrix2D Addm = new DenseDoubleMatrix2D(A);
//                        DenseDoubleMatrix2D mu0M = (DenseDoubleMatrix2D) a.subMatrix(Addm, 0, points.size() - 2, 0,
//                              points.size() - 2);
//
//                        double mu0 = a.det(mu0M);
//                        if (mu0 != 0)
//                        {
//                           for (int j = 0; j < m_dim - k; j += 1)
//                           {
//                              double bj = 0;
//                              double[] hj = new double[m_dim];
//                              for (int n = 0; n < m_dim; n += 1)
//                              {
//                                 if (n < k)
//                                 {
//                                    DenseDoubleMatrix2D mujn = (DenseDoubleMatrix2D) mu0M.copy(); // irgendwas aus mu0M
//                                    for (int i = 0; i < k; i += 1)
//                                    {
//                                       mujn.setQuick(i, n, A[i][k + j]);
//                                    }
//
//                                    hj[n] = a.det(mujn);
//                                 }
//                                 else
//                                    if (n == k + j)
//                                    {
//                                       hj[n] = -mu0;
//                                    }
//                                    else
//                                    {
//                                       hj[n] = 0;
//                                    }
//
//                                 bj = bj + 0.5 * hj[n] * (points.get(1)[n] + points.get(0)[n]);
//                              }
//
//                              A[k + j] = hj;
//                              B[k + j] = bj;
//                           }
//                           Addm = new DenseDoubleMatrix2D(A);
//                           if (a.det(Addm) != 0)
//                           {
//                              DenseDoubleMatrix1D Bddm = new DenseDoubleMatrix1D(B);
//                              cntr = a.mult(a.inverse(Addm), Bddm).toArray();
//                           }
//                           else
//                           {
//                              return null;
//                           }
//
//                        }
//                        else
//                        {
//                           return null;
//                           // throw new ArithmeticException("Not possible to create a Sphere based on given
//                           // coordinates.");
//                        }
//                     }
//                     else
                     {
                        double[][] mu0M = Base.getSubMatrix(A, 0, points.size() - 2, 0, points.size() - 2);
                        double mu0 = Base.determinant(mu0M);
                        if (mu0 != 0)
                        {
                           for (int j = 0; j < m_dim - k; j += 1)
                           {
                              double bj = 0;
                              double[] hj = new double[m_dim];
                              for (int n = 0; n < m_dim; n += 1)
                              {
                                 if (n < k)
                                 {
                                    // DenseDoubleMatrix2D mujn = (DenseDoubleMatrix2D) mu0M.copy(); // irgendwas aus
                                    // mu0M
                                    double[][] mujn = (double[][]) mu0M.clone();
                                    for (int i = 0; i < k; i += 1)
                                    {
                                       // mujn.setQuick(i, n, A[i][k + j]);
                                       mujn[i][n] = A[i][k + j];
                                    }

                                    hj[n] = Base.determinant(mujn);// a.det(mujn);
                                 }
                                 else
                                    if (n == k + j)
                                    {
                                       hj[n] = -mu0;
                                    }
                                    else
                                    {
                                       hj[n] = 0;
                                    }

                                 bj = bj + 0.5 * hj[n] * (points.get(1)[n] + points.get(0)[n]);
                              }

                              A[k + j] = hj;
                              B[k + j] = bj;
                           }
                           
                           if (m_dim == 2)
                              cntr = Base.solveExplicit2x2(A, B);
                           else
                              if (m_dim == 3)
                                 cntr = Base.solveExplicit3x3(A, B);
                              else
                                 cntr = Base.solve_Cramer(A, B, GMap.TOLERANCE_EPS);

                           if (cntr != null)
                           {
                              rds = Base.distance(cntr, points.get(0));
                           }
                           else
                              return null;
                        }
                        else
                           return null;
                     }

                  }
                  else
                  {
                     // TODO: replace this by a new implementation
                     // SphereN ns = Base.getCircumscribedNBall(points);
                     // cntr = ns.getCenterAsArray();
                     // rds = ns.getRadius();
                     return Base.getCircumscribedNBallByArr((ArrayList<double[]>) points);
                     // throw new ArithmeticException(
                     // "Not possible to create a Sphere based on given coordinates - to many coordinate points, not
                     // implemented.");
                  }
         return new Sphere(cntr, rds);
      }

      return null;
   }
}
