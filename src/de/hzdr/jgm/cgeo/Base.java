package de.hzdr.jgm.cgeo;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.IdentityHashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import de.hzdr.jgm.cgeo.gmap.GMap;
import de.hzdr.jgm.cgeo.gmap.Sphere;

/**
 * Base class for 2D/3D functionality. Encapsules fundamental functionality for 2D/3D computational geometry
 * 
 * Additionally: Provides some Matrix/LES functionality Provides some String 2 BaseTypeArray functionality for parsing.
 *
 * Note: GMap version, containing only functionality related to cgeo.gmap
 *
 * @author P. Menzel, Arbeitsgruppe Modellierung und Bewertung, HIF, HZDR
 * @version 0.0alpha
 */
public class Base
{

   /**
    * tolerance epsilon ... taken from the Lin.-Algebra-factory-member of Mat-class
    */
   static final double TOLERANCE_EPS = GMap.getInternalTolerance();// 0.000000001;// Math.pow(1,-9);
   
   // basic computational geometry

   public static boolean isEqualCoords(double[] coord1, double[] coord2)
   {
      return isEqualCoords(coord1, coord2, TOLERANCE_EPS);
   }

   public static boolean isEqualCoords(double[] coord1, double[] coord2, double tolerance)
   {
      boolean ret = false;

      if (coord1.length == coord2.length)
      {
         boolean is = true;
         int i = 0;

         while (is && i < coord1.length)
         {
            is = Math.abs(coord1[i] - coord2[i]) < tolerance;
            i++;
         }

         ret = is;
      }
      return ret;
   }

   /**
    * caluclates distance between to coordinates
    * 
    * @param pos1
    *           1st coordinate
    * @param pos2
    *           2nd coordinate
    * @return simple distance
    */
   public static double distance(double[] pos1, double[] pos2)
   {
      double[] diff = new double[pos1.length];
      for (int i = 0; i < pos1.length; i++)
         diff[i] = pos2[i] - pos1[i];

      return Base.norm(diff);
   }

   /**
    * calculates squared distance between to coordinates
    * 
    * @param pos1
    *           1st coordinate
    * @param pos2
    *           2nd coordinate
    * @return squared distance
    */
   public static double distanceSq(double[] pos1, double[] pos2)
   {
      double sqdiff = 0;
      for (int i = 0; i < pos1.length; i++)
         sqdiff = sqdiff + (pos2[i] - pos1[i]) * (pos2[i] - pos1[i]);

      return sqdiff;
   }

   /**
    * Calculates the length of a 2D vector.
    *
    * @param vec
    *           vector data as double array
    * @return length of the vector
    */
   public static double norm2(double[] vec)
   {
      try
      {
         return Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1]);
      }
      catch (ArrayIndexOutOfBoundsException e)
      {
         System.out.println("Warning: cgeo.Base.norm2(): " + e);
         return 0;
      }
   }

   /**
    * Calculates the length of a 3D vector.
    *
    * @param vec
    *           vector data as double array
    * @return length of the vector
    */
   public static double norm3(double[] vec)
   {
      return Math.sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
   }

   /**
    * Calculates the length of a n-dimensional vector.
    * 
    * Note: use norm2() or norm3 for n=2 or n=3
    *
    * @param vec
    *           vector data as double array
    * @return length of the vector
    */
   public static double norm(double[] vec)
   {
      if (vec.length == 2)
      {
         return norm2(vec);
      }
      else
         if (vec.length == 3)
         {
            return norm3(vec);
         }
         else
         {
            double ret = 0;
            for (int i = 0; i < vec.length; i += 1)
            {
               ret = ret + vec[i] * vec[i];
            }
            return Math.sqrt(ret);
         }
   }

   /**
    * Calculates the dot product of two 2D vectors.
    *
    * @param vec1
    *           first vector as double array
    * @param vec2
    *           second vector as double array
    * @return dot product of vec1 and vec2
    */
   public static double dot2(double[] vec1, double[] vec2)
   {
      try
      {
         return vec1[0] * vec2[0] + vec1[1] * vec2[1];
      }
      catch (ArrayIndexOutOfBoundsException e)
      {
         System.out.println("Warning: cgeo.Base.dot2(): " + e);
         return 0;
      }
   }

   /**
    * Calculates the dot product of two 3D vectors.
    *
    * @param vec1
    *           first vector as double array
    * @param vec2
    *           second vector as double array
    * @return dot product of vec1 and vec2
    */
   public static double dot3(double[] vec1, double[] vec2)
   {
      return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2];
   }

   /**
    * Calculates the standard dot product of two n-dimensional vectors.
    * 
    * Hint: use dot3 for 3d-vectors
    *
    * @param vec1
    *           first vector (double array)
    * @param vec2
    *           second vector (double array)
    * @return dot product of vec1 and vec2
    */
   public static double dot(double[] vec1, double[] vec2)
   {
      double ret = 0;
      if (vec1.length == vec2.length)
      {
         if (vec1.length == 2)
         {
            ret = dot2(vec1, vec2);
         }
         else
            if (vec1.length == 3)
            {
               ret = dot3(vec1, vec2);
            }
            else
            {
               for (int i = 0; i < vec1.length; i += 1)
               {
                  ret = ret + vec1[i] * vec2[i];
               }
            }
      }
      else
         throw new IllegalArgumentException("dimension missmatch");
      return ret;
   }

   /**
    * Calculates the cross product of two 3D vectors.
    *
    * @param vec1
    *           first vector as double array
    * @param vec2
    *           second vector as double array
    * @return cross product of vec1 and vec2 as double array
    */
   public static double[] cross3(double[] vec1, double[] vec2)
   {
      double[] cp =
      { 0, 0, 0 };
      try
      {
         cp[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1];
         cp[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2];
         cp[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0];
      }
      catch (ArrayIndexOutOfBoundsException e)
      {
         System.out.println("Warning: cgeo.Base.cross3(): " + e);
      }
      return cp;
   }

   // geometric operations

   /**
    * cuts 2 edges represented as start and end vertex
    * 
    * Note uses {@link Base#cutLineLine(double[], double[], double[], double[], double) cutLineLine}
    * 
    * @param e1s start point, first edge
    * @param e1e end point, first edge
    * @param e2s start point, second edge
    * @param e2e end point, second edge
    * @param tolerance
    * @return null for no cut or intersection point
    */
   public static double[] cutEdgeEdge(double[] e1s, double[] e1e, double[] e2s, double[] e2e, double tolerance)
   {
      int dim = e1s.length;
      double[] edge1 = new double[dim];
      double[] edge2 = new double[dim];
      
      for(int i = 0;i < dim; i++)
      {
         edge1[i] = e1e[i]-e1s[i];
         edge2[i] = e2e[i]-e2s[i];
      }
      
      double[] lambda = cutLineLine(e1s,edge1,e2s,edge2,tolerance);
      
      if(lambda[0] != Double.NaN && 
            lambda[0]>-tolerance && lambda[0]<(1+tolerance) &&
            lambda[1]>-tolerance && lambda[1]<(1+tolerance))
      {
         double[] interP = new double[dim];
         
         for(int i = 0;i < dim; i++)
         {
            interP[i] = e1s[i] + lambda[0]*edge1[i];
         }     
         return interP;
      }
      else
         return null;
   }

   /**
    * calculates cut between 2 lines , represented by a point (v) and a direction (dir)
    * 
    * @param v1
    *           first point
    * @param dir1
    *           first direction
    * @param v2
    *           second point
    * @param dir2
    *           second direction
    * @return 2 doubles containing the lambda values lamda == Nan for not cut
    */

   public static double[] cutLineLine(double[] v1, double[] dir1, double[] v2, double[] dir2, double tolerance)
   {
      double[] lambda =
      { Double.NaN, Double.NaN };

      
      
      double[] d = new double[v1.length];
      for (int i = 0; i < v1.length; i++)
         d[i] = v2[i] - v1[i];

      int x = 0;
      int y = 1;
      double denomimator = dir1[y] * dir2[x] - dir1[x] * dir2[y];
      
      if (Math.abs(denomimator) < tolerance && v1.length == 3)
      {
         if (Math.abs(dir1[x]) < tolerance && Math.abs(dir2[x]) < tolerance)
         {
            x = 2;
         }
         else
            if (Math.abs(dir1[y]) < tolerance && Math.abs(dir2[y]) < tolerance)
            {
               y = 2;
            }
            else
               x = 2;
         denomimator = dir1[y] * dir2[x] - dir1[x] * dir2[y];
      }
      else
      {
         // search adequate x/y for dim. > 3
      }
      
      if (Math.abs(denomimator) > tolerance)
      {
         if (dir2[x] != 0)
         {
            lambda[0] = (d[y] * dir2[x] - d[x] * dir2[y]) / denomimator;
            lambda[1] = (d[x] - lambda[0] * dir1[x]) / -dir2[x];
         }
         else
            if (dir1[x] != 0)
            {
               lambda[1] = (d[y] * dir1[x] - d[x] * dir1[y]) / denomimator;
               lambda[0] = (d[x] + lambda[1] * dir2[x]) / dir1[x];
            }

         if (lambda[0] != Double.NaN && lambda[1] != Double.NaN)
         {
            
            for (int i = 0; i < v1.length; i++)
            {
               if (i != x && i != y)
               {
                  if (Math.abs(v1[i] + lambda[0] * (dir1[i]) - v2[i]
                        - lambda[1] * (dir2[i])) > tolerance)
                  {
                     lambda[0] = Double.NaN;
                     lambda[1] = Double.NaN;
                     break;
                  }
               }
            }
         }
      }

      return lambda;
   }

//   /**
//    * cuts 2 edges represented as start and end vertex
//    * 
//    * Note uses {@link Base#cutLineLine_C(double[], double[], double[], double[]) cutLineLine_C}}
//    * 
//    * @param e1s start point, first edge
//    * @param e1e end point, first edge
//    * @param e2s start point, second edge
//    * @param e2e end point, second edge
//    * @param tolerance
//    * @return null for no cut or intersection point
//    */
//   public static double[] cutEdgeEdge_C(double[] e1s, double[] e1e, double[] e2s, double[] e2e, double tolerance)
//   {
//      int dim = e1s.length;
//      double[] edge1 = new double[dim];
//      double[] edge2 = new double[dim];
//      
//      for(int i = 0;i < dim; i++)
//      {
//         edge1[i] = e1e[i]-e1s[i];
//         edge2[i] = e2e[i]-e2s[i];
//      }
//      
//      double[] lambda = cutLineLine_C(e1s,edge1,e2s,edge2);
//      
//      
//      if(lambda[0] != Double.NaN && 
//            lambda[0]>-tolerance && lambda[0]<(1+tolerance) &&
//            lambda[1]>-tolerance && lambda[1]<(1+tolerance))
//      {
//         double[] interP = new double[dim];
//         
//         for(int i = 0;i < dim; i++)
//         {
//            interP[i] = e1s[i] + lambda[0]*edge1[i];
//         }     
//         return interP;
//      }
//      else
//         return null;
//   }
   
//   /**
//    * calculates cut between 2 lines , represented by a point (v) and a direction (dir)
//    * 
//    * Note: solve lambda1*dir1-lambda2*dir2 = v2-v1 with Colt
//    * 
//    * @param v1
//    *           first point
//    * @param dir1
//    *           first direction
//    * @param v2
//    *           second point
//    * @param dir2
//    *           second direction
//    * @return 2 doubles containing the lambda values lamda == Nan for not cut
//    */
//   public static double[] cutLineLine_C(double[] v1, double[] dir1, double[] v2, double[] dir2)
//   {
//
//      double[][] aArr = new double[dir1.length][2];
//      double[] bArr = new double[v1.length];
//
//      for (int i = 0; i < v1.length; i++)
//      {
//         aArr[i][0] = dir1[i];
//         aArr[i][1] = dir2[i];
//
//         bArr[i] = v2[i] - v1[i];
//      }
//
//      DenseDoubleMatrix2D A = new DenseDoubleMatrix2D(aArr);
//      DenseDoubleMatrix1D b = new DenseDoubleMatrix1D(bArr);
//      Algebra a = new Algebra();
//
//      double[] lambda = (a.mult(a.inverse(A), b)).toArray();
//      lambda[1] = -1 * lambda[1];
//
//      if (lambda[0] != Double.NaN && lambda[1] != Double.NaN)
//      {
//         for (int i = 0; i < v1.length; i++)
//         {
//            if (Math.abs(v1[i] + lambda[0] * (dir1[i]) - v2[i]
//                  - lambda[1] * (dir2[i])) > a.property().tolerance())
//            {
//               lambda[0] = Double.NaN;
//               lambda[1] = Double.NaN;
//               break;
//            }
//         }
//      }
//
//      return lambda;
//   }

   /**
    * Calculates the cut between a plane and an edges (or line) in 3D.
    * 
    * Plane: dot(plane.normal, Line: (X-planer.position) =0; ep1+lambda*(ep2-ep1)=X
    *
    * @param planeNormal
    *           normal of plane in 3D as double array
    * @param planePos
    *           point on plane in 3D as double array
    * @param ep1
    *           first point (as double array) of edge
    * @param ep2
    *           second point (as double array) of edge
    * @return Double lambda value. Both edge cuts, when lambda value is 0<=lambda<=1, for other lambda values the
    *         cutting position lies outside the edge on the defining line. When line and plane don't have a cutting
    *         point or infinity points, lambda = NaN.
    */
   public static double cutPlaneEdgeV2(double[] planeNormal, double[] planePos, double[] ep1, double[] ep2)
   {
      double lambda = Double.NaN;

      double d = (planeNormal[0] * (ep2[0] - ep1[0]) +
            planeNormal[1] * (ep2[1] - ep1[1]) +
            planeNormal[2] * (ep2[2] - ep1[2]));
      if (d != 0)
      {
         lambda = (-planeNormal[0] * ep1[0] + planeNormal[0] * planePos[0]
               - planeNormal[1] * ep1[1] + planeNormal[1] * planePos[1]
               - planeNormal[2] * ep1[2] + planeNormal[2] * planePos[2]) / d;
      }
      return lambda;
   }

   /**
    * caluclates the intersecting edge for two flat non-parallel 3D polygons
    * 
    * @param poly1
    *           ordered seq. of vertex coordiantes for first polygon
    * @param poly2
    *           ordered seq. of vertex coordinates for second polygon
    * @return two edge points
    */
   public static double[][] cutFlatPolyPoly(double[][] poly1, double[][] poly2, double tolerance)
   {
      // workflow:
      // 1) calculates intersecting line between 2 hyperplanes
      // 2) intersect line with poly edges
      int dim = poly1[0].length;
      if (dim == 3)
      {
         double[][] edgePoints = new double[2][dim];
         double[] e11 = new double[dim];
         double[] e12 = new double[dim];
         double[] e21 = new double[dim];
         double[] e22 = new double[dim];

         for (int i = 0; i < dim; i++)
         {
            e11[i] = poly1[1][i] - poly1[0][i];
            e12[i] = poly1[poly1.length - 1][i] - poly1[0][i];

            e21[i] = poly2[1][i] - poly2[0][i];
            e22[i] = poly2[poly2.length - 1][i] - poly2[0][i];
         }        
         
         double[] normal1 = cross3(e11, e12);
         double[] normal2 = cross3(e21, e22);
         
         
         if (Math.abs(dot3(normal1, normal2)) - 1 < tolerance)
         {
            // Intersection of two planes:
            // https://de.wikipedia.org/wiki/Schnittgerade#Schnitt_zweier_Ebenen_in_Normalenform
            // there may be better solutions

            double[] w = cross3(normal1, normal2);
            
            double d1 = dot3(normal1, poly1[0]);
            double d2 = dot3(normal2, poly2[0]);

            double[] q = new double[dim];

            double n1Sq = dot3(normal1, normal1);
            double n2Sq = dot3(normal2, normal2);
            double n1n2 = dot3(normal1, normal2);

            double denominator = (n1Sq * n2Sq - n1n2 * n1n2);

            double coeff1 = (d1 * n2Sq - d2 * n1n2) / denominator;
            double coeff2 = (d2 * n1Sq - d1 * n1n2) / denominator;

            for (int i = 0; i < dim; i++)
            {
               q[i] = coeff1 * normal1[i] + coeff2 * normal2[i];
            }
            
            // intersect ray q+w with all edges
            // collect all ray lambda for 0<=edge-lambda<=1
            
            double[] l1 = new double[poly1.length];
            int cN1 = 0;
            
            Map<Double,Integer> lambdaValues = new IdentityHashMap<Double,Integer>();
            
            for (int i = 0; i < poly1.length; i++)
            {
               double[] p = poly1[i];
               double[] edge = new double[dim];
               
               int nextId = i+1;
               if(i == poly1.length-1)
                  nextId = 0;
               
               for (int j = 0; j < dim; j++)
               {
                  edge[j] = poly1[nextId][j]-p[j];
               }
               
               double[] lambda = cutLineLine(q,w,p,edge,tolerance);
               
               if(lambda[1]>-tolerance && lambda[1]<1+tolerance)
               {
                  //todo save lambda[0]
//                  l1[cN1] = lambda[0];
//                  cN1++;
                  System.out.println("poly 1, edge "+i+": Lambda values: "+Arrays.toString(lambda));
                  lambdaValues.put(lambda[0], 0);
               }
            }

            double[] l2 = new double[poly1.length];
            int cN2 = 0;     
            
            for (int i = 0; i < poly2.length; i++)
            {
               double[] p = poly2[i];
               double[] edge = new double[dim];
               
               int nextId = i+1;
               if(i == poly2.length-1)
                  nextId = 0;
               
               for (int j = 0; j < dim; j++)
               {
                  edge[j] = poly2[nextId][j]-p[j];
               }
               double[] lambda = cutLineLine(q,w,p,edge,tolerance);
               if(lambda[1]>-tolerance && lambda[1]<1+tolerance)
               {
                  //todo save lambda[0]
//                  l2[cN2] = lambda[0];
//                  cN2++;
                  System.out.println("poly 2, edge "+i+": Lambda values: "+Arrays.toString(lambda));
                  lambdaValues.put(lambda[0], 1);
               }
            }
            
            //for convex polygons there should be zero or 2 lambda per polygon
            //sort all lambdas
            //

            System.out.println("Lambda values: "+lambdaValues);
            
            if(lambdaValues.size()==4)
            {                             
//               double[] lambda = new double[4];
               
               //sort lambda here
               
               ArrayList<Double> lambdas2sort = new ArrayList<Double>(lambdaValues.keySet());
               
               Collections.sort(lambdas2sort);
               
               if(lambdaValues.get(lambdas2sort.get(0)) != lambdaValues.get(lambdas2sort.get(1)))
               {
                  for (int j = 0; j < dim; j++)
                  {
                     edgePoints[0][j] = q[j]+((double)lambdas2sort.get(1))*w[j];
                     edgePoints[1][j] = q[j]+((double)lambdas2sort.get(2))*w[j];
                  }
               
                  return edgePoints;
               }
               else
                  return null;
            }
            else
               return null;
         }
         else
         {
            return null;
         }
      }
      else
      {
         return null;
      }

   }

   /**
    * calculates the minimum circumscribed N-ball for a Set of n-dimensional points
    * 
    * @param points
    * @return
    */
   @SuppressWarnings("unchecked")
   public static Sphere getCircumscribedNBallByArr(ArrayList<double[]> points)
   {
      ArrayList<double[]> C = (ArrayList<double[]>) points.clone();
      ArrayList<double[]> B = new ArrayList<double[]>();

      return minCircSphByArr(C, B, (ArrayList<double[]>) C.clone());
   }  
   /**
    * calculates recursively the minimum circumscribed N-ball for a Set of N-dimensional points
    * 
    * implementation of Xu et al. (2003): Solution Methodologies for the Smallest Enclosing Circle Problem
    * 
    * Welzl (1991): Smallest enclosing disks (balls and ellipsoids)
    * 
    * @param points
    * @return
    */
   private static Sphere minCircSphByArr(ArrayList<double[]> C, ArrayList<double[]> B, ArrayList<double[]> orig)
   {
      Sphere D = new Sphere();

      if (C.size() == 0)
      {
         if (B.size() == 1)
         {
            D = new Sphere(B.get(0));
         }
         else
            if (B.size() == 2)
            {
               D = new Sphere(B.get(0), B.get(1));
            }
            else
               if (B.size() > 2)
               {
                  D.defineBy(B);
               }
      }
      else
      {
         Random rndGen = new Random();
         int index = rndGen.nextInt(C.size());
         double[] p = C.get(index);
         C.remove(index);
         D = minCircSphByArr(C, B, orig);
         if (D != null && !D.isIn(p))
         {
            B.add(p);
            D = minCircSphByArr(C, B, orig);
         }
      }

      return D;
   }   

   public static double[] getFacetNormalNormalized(double[][] points)
   {
      // computing the normalized normal
      double[] ps =
      { points[1][0] - points[0][0], points[1][1] - points[0][1], points[1][2] - points[0][2] };
      double[] pe =
      { points[points.length - 1][0] - points[0][0], points[points.length - 1][1] - points[0][1],
            points[points.length - 1][2] - points[0][2] };
      double[] normal = Base.cross3(ps, pe);

      double nNorm = Base.norm3(normal);
      normal[0] = normal[0] / nNorm;
      normal[1] = normal[1] / nNorm;
      normal[2] = normal[2] / nNorm;

      return normal;
   }

   /**
    * Area of a planar polygon
    * 
    * citation: R. N. Goldman. IV.1 - AREA OF PLANAR POLYGONS AND VOLUME OF POLYHEDRA. In Graphics Gems II
    * (ed. J. ARVO), S. 170 – 171. Morgan Kaufmann, San Diego (1991). ISBN 978-0-08-050754-5.
    * doi:http://dx.doi.org/10.1016/B978-0-08-050754-5.50043-8.
    * 
    * 
    * @param points
    *           set of points defining the polygons border
    * @return area
    */
   public static double polygonArea3(double[][] points)
   {

      double[] normal = Base.getFacetNormalNormalized(points);

      double[] area = new double[3];

      for (int i = 0; i < points.length; i++)
      {
         double[] cp;
         if (i < points.length - 1)
            cp = Base.cross3(points[i], points[i + 1]);
         else
            cp = Base.cross3(points[i], points[0]);
         area[0] = area[0] + cp[0];
         area[1] = area[1] + cp[1];
         area[2] = area[2] + cp[2];
      }

      return 0.5 * Math.abs(Base.dot(normal, area));
   }

   /**
    * Area of a planar polygon
    * 
    * citation: R. N. Goldman. IV.1 - AREA OF PLANAR POLYGONS AND VOLUME OF POLYHEDRA. In Graphics Gems II
    * (ed. J. ARVO), S. 170 – 171. Morgan Kaufmann, San Diego (1991). ISBN 978-0-08-050754-5.
    * doi:http://dx.doi.org/10.1016/B978-0-08-050754-5.50043-8.
    * 
    * 
    * @param points
    *           set of points defining the polygons border
    * @param normal
    *           facet's normalized normal
    * @return area
    */
   public static double polygonArea3(double[][] points, double[] normal)
   {

      double[] area = new double[3];

      for (int i = 0; i < points.length; i++)
      {
         double[] cp;
         if (i < points.length - 1)
            cp = Base.cross3(points[i], points[i + 1]);
         else
            cp = Base.cross3(points[i], points[0]);
         area[0] = area[0] + cp[0];
         area[1] = area[1] + cp[1];
         area[2] = area[2] + cp[2];
      }

      return 0.5 * Math.abs(Base.dot(normal, area));
   }

   /**
    * volume of a polyhedron with planar facets
    * 
    * citation: R. N. Goldman. IV.1 - AREA OF PLANAR POLYGONS AND VOLUME OF POLYHEDRA. In Graphics Gems II
    * (ed. J. ARVO), S. 170 – 171. Morgan Kaufmann, San Diego (1991). ISBN 978-0-08-050754-5.
    * doi:http://dx.doi.org/10.1016/B978-0-08-050754-5.50043-8.
    * 
    * 
    * @param facets
    *           set of facet defininitions
    * @return area
    */
   public static double polyhedronVolume(double[][][] facets)
   {
      double volume = 0;

      for (int f = 0; f < facets.length; f++)
      {
         double[] n = Base.getFacetNormalNormalized(facets[f]);
         double facetArea = Base.polygonArea3(facets[f], n);

         volume = volume + (Base.dot3(facets[f][0], n)) * facetArea;
      }
      return (1 / 3) * Math.abs(volume);
   }

   /**
    * surface area of a polyhedron with planar facets
    * 
    * citation: R. N. Goldman. IV.1 - AREA OF PLANAR POLYGONS AND VOLUME OF POLYHEDRA. In Graphics Gems II
    * (ed. J. ARVO), S. 170 – 171. Morgan Kaufmann, San Diego (1991). ISBN 978-0-08-050754-5.
    * doi:http://dx.doi.org/10.1016/B978-0-08-050754-5.50043-8.
    * 
    * 
    * @param facets
    *           set of facet defininitions
    * @return area
    */
   public static double polyhedronSurfaceArea(double[][][] facets)
   {
      double sarea = 0;

      for (int f = 0; f < facets.length; f++)
      {
         sarea = sarea + Base.polygonArea3(facets[f]);
      }
      return sarea;
   }

   /**
    * Area of a planar polygon
    * 
    * citation: R. N. Goldman. IV.1 - AREA OF PLANAR POLYGONS AND VOLUME OF POLYHEDRA. In Graphics Gems II
    * (ed. J. ARVO), S. 170 – 171. Morgan Kaufmann, San Diego (1991). ISBN 978-0-08-050754-5.
    * doi:http://dx.doi.org/10.1016/B978-0-08-050754-5.50043-8.
    * 
    * 
    * @param points
    *           set of points defining the polygons border
    * @return area
    */
   public static double polygonArea2(double[][] points)
   {

      double[] area = new double[3];

      for (int i = 0; i < points.length; i++)
      {
         double[] cp;
         if (i < points.length - 1)
            cp = Base.cross3(new double[]
            { points[i][0], points[i][1], 0 }, new double[]
            { points[i + 1][0], points[i + 1][1], 0 });
         else
            cp = Base.cross3(new double[]
            { points[i][0], points[i][1], 0 }, new double[]
            { points[0][0], points[0][1], 0 });
         area[0] = area[0] + cp[0];
         area[1] = area[1] + cp[1];
         area[2] = area[2] + cp[2];
      }

      return 0.5 * Base.norm3(area);
   }

   // basic Matrix/LES
   /**
    * Calculates the determinant of 2x2 matrix.
    * 
    * @param A
    *           Matrix as double[2][2] array
    * 
    * @return determinant of A
    */
   public static double determinant2(double[][] A)
   {
      if (A[0].length >= 2 && A.length >= 2)
         return (A[0][0] * A[1][1] - A[1][0] * A[0][1]);
      else
         return 0;
   }

   /**
    * Calculates the determinant of 3x3 matrix.
    * 
    * (based on MatrixTools.h (C++) by P. Menzel)
    * 
    * @param A
    *           Matrix as double[3][3] array
    * 
    * @return determinant of A
    */
   public static double determinant3(double[][] A)
   {
      if (A[0].length >= 3 && A.length >= 3)
         return (A[2][2] * A[1][1] - A[2][1] * A[1][2]) * A[0][0] +
               (-A[2][2] * A[0][1] + A[2][1] * A[0][2]) * A[1][0] +
               (-A[0][2] * A[1][1] + A[0][1] * A[1][2]) * A[2][0];
      else
         return 0;
   }

   /**
    * Calculates the determinant of nxn matrix.
    * 
    * Hint: for n=2 or n=3 use determinant2 or determinant 3
    * 
    * @param A
    *           Matrix as double[n][n] array
    * @return determinant of A
    */
   public static double determinant(double[][] A)
   {
      double det = 0;
      if (A.length == 2)
      {
         det = determinant2(A);
      }
      else
         if (A.length == 3)
         {
            det = determinant3(A);
         }
         else
         {
            // TODO: develop after Laplace Development Rule using sub-determinants
            int r = 0; // development after first row, may be adapted later

            for (int c = 0; c < A[r].length; c += 1)
            {
               if (A[r][c] != 0)
               {
                  det = det + (A[r][c] * Math.pow(-1, r + c) * determinant(getSubMatrix(A, r, c)));
               }
            }
         }
      return det;
   }

//   /**
//    * Calculates the determinant of a matrix by using Colt.linAlg.Algebra functionalities
//    * 
//    * Note: this is only for example, how we may use Colt capsuled
//    * 
//    * @param A
//    *           Matrix as double[n][n] array
//    * @return determinant of A
//    */
//   public static double determinant_C(double[][] A)
//   {
//      Algebra a = new Algebra();
//      return a.det(new DenseDoubleMatrix2D(A));
//   }

   /**
    * Creates a submatrix (n-1)x(m-1) for matrix A without row i and column j
    * 
    * @param A
    *           Matrix as double[n][m] array
    * @param i
    *           row
    * @param j
    *           column
    * @return submatrix as double[n-1][m-1] array
    */
   public static double[][] getSubMatrix(double[][] A, int i, int j)
   {
      if (A.length >= 3 && A[0].length >= 3)
      {
         double[][] ret = new double[A.length - 1][A[0].length];
         int nr = 0;
         for (int r = 0; r < A.length; r += 1)
         {
            if (r != i)
            {
               if (j > 0)
                  System.arraycopy(A[r], 0, ret[nr], 0, j);
               if (j < (A[0].length - 1))
                  System.arraycopy(A[r], j + 1, ret[nr], j + 1, A[0].length - j - 1);
               nr += 1;
            }
         }

         return ret;
      }
      else
         return null;
   }

   public static double[][] getSubMatrix(double[][] A, int i1, int i2,int j1, int j2)
   {
      double[][] ret = new double[(i2-i1)+1][(j2-j1)+1];
      
      for(int i = i1; i <=i2;i++)
      {
         for(int j = j1; j <=j2;j++)
         {
            ret[i][j] = A[i][j];
         }         
      }
      return ret;
   }

   /**
    * Solves the linear equation system Ax = b by using an explicit solver
    * 
    * 
    * @param A
    *           Matrix as double[2][2] array
    * @param b
    *           right side of the system as double[2] array
    * 
    * @return x as double[2] array or null if determinant(A)=0
    */
   public static double[] solveExplicit2x2(double[][] A, double[] b)
   {
      double det = determinant2(A);
      double[] x = new double[2];
      if(det != 0)
      {
         double nen = 1.0 / det;
         x[0] = (A[1][1]*b[0]-A[0][1]*b[1])*nen;
         x[1] = (A[0][0]*b[1]-A[1][0]*b[0])*nen;
         return x;
      }
      else
         return null;
   }
   
   /**
    * Solves the linear equation system Ax = b by using an explicit solver
    * 
    * (based on MatrixTools.h (C++) by P. Menzel)
    * 
    * @param A
    *           Matrix as double[3][3] array
    * @param b
    *           right side of the system as double[3] array
    * 
    * @return x as double[3] array or null if determinant(A)=0
    */
   public static double[] solveExplicit3x3(double[][] A, double[] b)
   {
      double det = determinant3(A);
      double[] x = new double[3];
      if (det != 0)
      {
         double nen = 1.0 / det;
         x[0] = (A[0][1] * A[1][2] * b[2] - A[0][1] * b[1] * A[2][2] +
               A[0][2] * A[2][1] * b[1] - A[0][2] * b[2] * A[1][1] +
               b[0] * A[2][2] * A[1][1] - b[0] * A[2][1] * A[1][2]) * nen;
         x[1] = -((A[1][2] * b[2] - b[1] * A[2][2]) * A[0][0] +
               (-A[0][2] * b[2] + b[0] * A[2][2]) * A[1][0] +
               (b[1] * A[0][2] - A[1][2] * b[0]) * A[2][0]) * nen;
         x[2] = ((-A[2][1] * b[1] + b[2] * A[1][1]) * A[0][0] +
               (-b[2] * A[0][1] + A[2][1] * b[0]) * A[1][0] +
               (-b[0] * A[1][1] + A[0][1] * b[1]) * A[2][0]) * nen;
         return x;
      }
      else
         return null;

   }

   /**
    * solves a n x x equation system with Cramer's rule, not efficient, but stable
    * @param A
    * @param b
    * @return x as double[n] array or null if determinant(A)=0
    */
   public static double[] solve_Cramer(double[][] A, double[] b,double tolerance)
   {
      double[] x = new double[b.length];
      
      double detA = determinant(A);
      
      if(Math.abs(detA)>tolerance)
      {
         
         for(int i = 0; i <x.length;i++)
         {
            //build Ai           
            double[][] Ai = new double[A.length][A[0].length];
            
            for(int v = 0; v < Ai.length;v++)
            {
               for(int u = 0; u < Ai[0].length;u++)
               {
                  if(u==i)
                  {
                     Ai[v][u] = b[v];
                  }
                  else
                  {
                     Ai[v][u] = A[v][u];
                  }
               }  
               double detAi = determinant(Ai);
               x[i] = detAi/detA;
            }
         }
         return x;
      }
      else
         return null;
      
   }
   
   /**
    * Creates all combinations of m indices out of n indices
    * 
    * Note: This method works, but seems to be very inefficient based on
    * http://www.java-forum.org/thema/alle-kombinationen-aus-arraylist-potenzmenge.133218/
    * 
    * @param numberOfElements
    *           number of indices
    * @param numberOfElementsPerKombuination
    *           number of indices per combination
    * 
    * @return ArryList of all combinations
    */
   public static ArrayList<List<Integer>> getAllIndexCombinations(int numberOfElements,
         int numberOfElementsPerKombuination)
   {
      return getAllIndexCombinations(numberOfElements, numberOfElementsPerKombuination, false);
   }

   /**
    * Creates all combinations of m indices out of n indices
    * 
    * Note: This method works, but seems to be very inefficient; based on
    * http://www.java-forum.org/thema/alle-kombinationen-aus-arraylist-potenzmenge.133218/
    * 
    * Note: special cases for combinations of 2,3,4 entries, allow higher number of entries Note: general case only
    * apply to arrays with few entriess (< 25)
    * 
    * @param numberOfElements
    *           number of indices
    * @param numberOfElementsPerKombuination
    *           number of indices per combination
    * @param doSort
    *           if true, perform a sort for combination length
    * 
    * @return ArryList of all combinations
    */
   public static ArrayList<List<Integer>> getAllIndexCombinations(int numberOfElements,
         int numberOfElementsPerKombuination, boolean doSort)
   {
      ArrayList<List<Integer>> kombinations = new ArrayList<List<Integer>>();
      int[] ids = new int[numberOfElements];
      for (int i = 0; i < numberOfElements; i += 1)
      {
         ids[i] = i;
      }

      if (numberOfElementsPerKombuination > numberOfElements)
         numberOfElementsPerKombuination = numberOfElements;

      if (numberOfElementsPerKombuination == 2)
      {
         for (int i = 0; i < numberOfElements; i += 1)
         {
            for (int j = i + 1; j < numberOfElements; j += 1)
            {
               List<Integer> subList = Arrays.asList(i, j);
               kombinations.add(subList);
            }
         }
      }
      else
         if (numberOfElementsPerKombuination == 3)
         {
            for (int i = 0; i < numberOfElements; i += 1)
            {
               for (int j = i + 1; j < numberOfElements; j += 1)
               {
                  for (int k = j + 1; k < numberOfElements; k += 1)
                  {
                     List<Integer> subList = Arrays.asList(i, j, k);
                     kombinations.add(subList);
                  }
               }
            }
         }
         else
            if (numberOfElementsPerKombuination == 4)
            {
               for (int i = 0; i < numberOfElements; i += 1)
               {
                  for (int j = i + 1; j < numberOfElements; j += 1)
                  {
                     for (int k = j + 1; k < numberOfElements; k += 1)
                     {
                        for (int l = k + 1; l < numberOfElements; l += 1)
                        {
                           List<Integer> subList = Arrays.asList(i, j, k, l);
                           kombinations.add(subList);
                        }
                     }
                  }
               }
            }
            else
            {

               for (long i = 0; i < 1 << numberOfElements; i++)
               {
                  List<Integer> subList = new ArrayList<Integer>();

                  for (int k = 0; k < numberOfElements; k++)
                  {

                     if ((i & (1 << k)) > 0)
                     {
                        subList.add(ids[k]);
                     }
                  }
                  if (numberOfElementsPerKombuination <= 0 || subList.size() == numberOfElementsPerKombuination)
                  {
                     kombinations.add(subList);
                  }
               }

               if (numberOfElementsPerKombuination <= 0 && doSort)
               {
                  // idea from http://javatricks.de/tricks/objekte-vergleichen-mit-comparable-und-comparator#header-1
                  Collections.sort(kombinations, new Comparator<List<Integer>>()
                  {
                     @Override
                     public int compare(List<Integer> first, List<Integer> second)
                     {
                        if (first.size() < second.size())
                           return -1;
                        else
                           if (first.size() == second.size())
                              return 0;
                           else
                              return 1;
                     }
                  });
               }
            }
      return kombinations;

   }
}
