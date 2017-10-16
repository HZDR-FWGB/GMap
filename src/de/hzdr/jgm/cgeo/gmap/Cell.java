package de.hzdr.jgm.cgeo.gmap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import de.hzdr.jgm.cgeo.Base;

/**
 * new realization for cell embedding
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class Cell implements CellObject
{

   private int      dimension;
   private Dart     cellDart;

   private Object   data   = null;

   /*
    * if dimension = 0, thsi is a vertex embedding
    */
   private double[] vertex = null;

   public Cell(int dim)
   {
      dimension = dim;
      vertex = new double[3];
   }

   public Cell(double... pointCoordinates)
   {
      dimension = 0;
      vertex = pointCoordinates;
   }

   public Dart getDart()
   {
      return cellDart;
   }

   public void setDart(Dart dart)
   {
      cellDart = dart;
   }

   /**
    * returns the reference to the underlying double array
    * @return reference to coordinate array
    */
   public double[] getVertex()
   {
      return this.vertex;
   }

   public void setVertex(double... pointCoordinates)
   {
      vertex = pointCoordinates;
   }

   public void setData(Object dataObj)
   {
      this.data = dataObj;
   }

   public Object getData()
   {
      return this.data;
   }

   // interface functionality, not overridden

   /**
    * Are these coordinates in or on the cell
    * 
    * @param pointCoordinates
    * @return
    */
   public boolean containsPoint(double... pointCoordinates)
   {
      return containsPoint(GMap.TOLERANCE_EPS, pointCoordinates);
   }

   /**
    * Are these coordinates in or on the cell
    * 
    * @param tolerance
    * @param pointCoordinates
    * @return
    */
   public boolean containsPoint(double tolerance, double... pointCoordinates)
   {
      boolean ret = false;
      if (dimension == 0)
      {
         // System.out.println(">>containsPoint point in point");
         return Base.isEqualCoords(pointCoordinates, vertex, tolerance);
      }
      else
         if (dimension == 1)
         {
            // System.out.println(">>containsPoint point in edge");
            // TODO: point on edge test
            Cell vertex1 = this.getDart().getVertex();
            Cell vertex2 = this.getDart().getInvolution(0).getVertex();

            if (vertex1.containsPoint(tolerance, pointCoordinates)
                  || vertex2.containsPoint(tolerance, pointCoordinates))
            {
               return true;
            }
            else
            {
               // Workflow:
               // test both edge vertices to contain the test point
               // calc. the edge vector and the vector between start-point and test coordinates
               // if scalar product is ~1 and edge is longer than test vector, then true
               double[] edgeVertex1 = vertex1.getVertex();
               double[] edgeVertex2 = vertex2.getVertex();

               double[] edge = new double[edgeVertex1.length];
               double[] test = new double[edgeVertex1.length];

               double n1 = 0;
               double n2 = 0;

               for (int i = 0; i < edgeVertex1.length; i++)
               {
                  edge[i] = edgeVertex2[i] - edgeVertex1[i];
                  test[i] = pointCoordinates[i] - edgeVertex1[i];

                  n1 = n1 + (edge[i] * edge[i]);
                  n2 = n2 + (test[i] * test[i]);
               }

               if (Math.abs(Base.dot(edge, test) - 1) < tolerance && (n2 - n1) < tolerance)
                  return true;
               else
                  return false;
            }
         }
         else
            if (dimension == 2)
            {
//               System.out.println(">>containsPoint point: "+Arrays.toString(pointCoordinates)+" in facet");
               // TODO: scanline test for edges
               // 0) test all vertices to contain the test point?
               // 1) define ray direction (testpoint - mean?)
               // 2) test each edge for cut with ray
               // 3) odd number of cuts -> cut!

//                CellList vertices = getVertices();
//                for(int i =0; i < vertices.size();i++)
//                   System.out.println(">>:"+Arrays.toString(vertices.get(i).vertex));
               // double[] center = vertices.getCenterOfMass();
               double[] center = getCenterOfMass();
               double[] ray = new double[center.length];
               if (Base.isEqualCoords(center, pointCoordinates))
               {
                  for (int i = 0; i < center.length; i++)
                  {
                     ray[i] = Math.random();
                  }
               }
               else
               {
                  for (int i = 0; i < center.length; i++)
                  {
                     ray[i] = center[i] - pointCoordinates[i];
                  }
               }
//               System.out.println(">>> ray:"+Arrays.toString(ray));

               Dart d = cellDart;
               Dart nd = d.getInvolution(0);
               int nCuts = 0; // number of edges cut
               while (nd != d.getInvolution(1))
               {
                  double[] v1 = nd.getVertexCoordinates();
                  double[] v2 = nd.getInvolution(0).getVertexCoordinates();

                  double[] edge = new double[v1.length];

                  for (int i = 0; i < v1.length; i++)
                  {
                     edge[i] = v2[i] - v1[i];
                  }
                  // System.out.println(">>>> edge: "+Arrays.toString(edge));

                  double lambda[] = Base.cutLineLine(pointCoordinates, ray, v1, edge, tolerance);
                  // System.out.println(">>>> lambdas: "+Arrays.toString(lambda));
                  // System.out.println((lambda[1] > tolerance)+"/"+(lambda[1] < (1 + tolerance)));
                  if (lambda[0] > 0 && lambda[1] > tolerance && lambda[1] < (1 + tolerance))
                     nCuts++;

                  nd = nd.getInvolution(1).getInvolution(0);

               }

               // test last edge
               double[] v1 = nd.getVertexCoordinates();
               double[] v2 = nd.getInvolution(0).getVertexCoordinates();

               double[] edge = new double[v1.length];

               for (int i = 0; i < v1.length; i++)
               {
                  edge[i] = v2[i] - v1[i];
               }
               // System.out.println(">>>> edge: "+Arrays.toString(edge));
               double lambda[] = Base.cutLineLine(pointCoordinates, ray, v1, edge, tolerance);
               // System.out.println(">>>> lambdas: "+Arrays.toString(lambda));
               // System.out.println((lambda[1] > tolerance)+"/"+(lambda[1] < (1 + tolerance)));
               if (lambda[0] > 0 && lambda[1] > tolerance && lambda[1] < (1 + tolerance))
                  nCuts++;

//               System.out.println(">>> number of cuts:"+nCuts);

               if (nCuts % 2 != 0)
               {
                  return true; // odd number of cuts
               }
               else
                  return false; // equal number of cuts

            }
            else
               if (dimension == 3)
               {
//                  System.out.println(">>containsPoint point in polyhedron");
                  // TODO: scanline test for facets
                  CellList vertices = getVertices();
                  double[] center = vertices.getCenterOfMass();

                  double[] ray = new double[center.length];
                  for (int i = 0; i < center.length; i++)
                  {
                     ray[i] = center[i] - pointCoordinates[i];
                  }

                  CellList facets = getBorderingCells(2);
                  int nCuts = 0; // number of edges cut
                  for (int f = 0; f < facets.size(); f++)
                  {
                     Cell facet = facets.get(f);
                     double[] p1 = facet.getDart().getVertexCoordinates();
                     double[] p2 = facet.getDart().getInvolution(0)
                           .getVertexCoordinates();
                     double[] p3 = facet.getDart().getInvolution(1)
                           .getInvolution(0)
                           .getVertexCoordinates();

                     double[] edge1 = new double[p1.length];
                     double[] edge2 = new double[p1.length];

                     for (int i = 0; i < p1.length; i++)
                     {
                        edge1[i] = p2[i] - p1[i];
                        edge2[i] = p3[i] - p1[i];
                     }

                     double[] facetNormal = Base.cross3(edge1, edge2);

                     double lambda = Base.cutPlaneEdgeV2(facetNormal, p1, pointCoordinates, ray);

                     if (lambda != Double.NaN)
                     {
                        double[] pointOnPlane = new double[p1.length];
                        for (int i = 0; i < p1.length; i++)
                        {
                           pointOnPlane[i] = pointCoordinates[i] + lambda * ray[i];
                        }

                        if (facet.containsPoint(tolerance, pointOnPlane))
                        {
                           nCuts++;
                        }
                     }
                  }
                  
//                  System.out.println("number of cuts: "+nCuts);
                  
                  if (nCuts % 2 != 0)
                  {
                     return true; // odd number of cuts
                  }
                  else
                     return false; // equal number of cuts

               }
               else
               {
                  System.out.println("Cell.containsPoint: not implemented for dimension = " + dimension + "!");
               }
      return ret;
   }

   /**
    * Are these coordinates in or on the cell
    * 
    * Note: version with external given center point
    * 
    * @param tolerance
    * @param center
    * @param pointCoordinates
    * @return
    */
   public boolean containsPoint(double tolerance, double[] center, double... pointCoordinates)
   {
      boolean ret = false;
      if (dimension == 0)
      {
         return Base.isEqualCoords(pointCoordinates, vertex, tolerance);
      }
      else
         if (dimension == 1)
         {
            // TODO: point on edge test
            Cell vertex1 = this.getDart().getVertex();
            Cell vertex2 = this.getDart().getInvolution(0).getVertex();

            if (vertex1.containsPoint(tolerance, pointCoordinates)
                  || vertex2.containsPoint(tolerance, pointCoordinates))
            {
               return true;
            }
            else
            {
               // Workflow:
               // test both edge vertices to contain the test point
               // calc. the edge vector and the vector between start-point and test coordinates
               // if scalar product is ~1 and edge is longer than test vector, then true
               double[] edgeVertex1 = vertex1.getVertex();
               double[] edgeVertex2 = vertex2.getVertex();

               double[] edge = new double[edgeVertex1.length];
               double[] test = new double[edgeVertex1.length];

               double n1 = 0;
               double n2 = 0;

               for (int i = 0; i < edgeVertex1.length; i++)
               {
                  edge[i] = edgeVertex2[i] - edgeVertex1[i];
                  test[i] = pointCoordinates[i] - edgeVertex1[i];

                  n1 = n1 + (edge[i] * edge[i]);
                  n2 = n2 + (test[i] * test[i]);
               }

               if (Math.abs(Base.dot(edge, test) - 1) < tolerance && (n2 - n1) < tolerance)
                  return true;
               else
                  return false;
            }
         }
         else
            if (dimension == 2)
            {
               // TODO: scanline test for edges
               // 0) test all vertices to contain the test point?
               // 1) define ray direction (testpoint - mean?)
               // 2) test each edge for cut with ray
               // 3) odd number of cuts -> cut!

//               CellList vertices = getVertices();
               // double[] center = vertices.getCenterOfMass();

               double[] ray = new double[center.length];
               for (int i = 0; i < center.length; i++)
               {
                  ray[i] = center[i] - pointCoordinates[i];
               }

               Dart d = cellDart;
               Dart nd = d.getInvolution(0);
               int nCuts = 0; // number of edges cut
               while (nd != d.getInvolution(1))
               {
                  double[] v1 = nd.getVertexCoordinates();
                  double[] v2 = nd.getInvolution(0).getVertexCoordinates();

                  double[] edge = new double[v1.length];

                  for (int i = 0; i < v1.length; i++)
                  {
                     edge[i] = v2[i] - v1[i];
                  }

                  // TODO: cutEdgeRay -> lambda for edge
                  // 0<=lambda<=1 cut++;

                  double lambda[] = Base.cutLineLine(pointCoordinates, ray, v1, edge, tolerance);

                  if (lambda[0] > 0 && lambda[1] > -tolerance && lambda[1] < (1 + tolerance))
                     nCuts++;

                  nd = nd.getInvolution(1).getInvolution(0);

               }

               if (nCuts % 2 != 0)
               {
                  return true; // odd number of cuts
               }
               else
                  return false; // equal number of cuts

            }
            else
               if (dimension == 3)
               {
                  // TODO: scanline test for facets
//                  CellList vertices = getVertices();
                  // double[] center = vertices.getCenterOfMass();

                  double[] ray = new double[center.length];
                  for (int i = 0; i < center.length; i++)
                  {
                     ray[i] = center[i] - pointCoordinates[i];
                  }

                  CellList facets = getBorderingCells(2);
                  int nCuts = 0; // number of edges cut
                  for (int f = 0; f < facets.size(); f++)
                  {
                     Cell facet = facets.get(f);
                     double[] p1 = facet.getDart().getVertexCoordinates();
                     double[] p2 = facet.getDart().getInvolution(0)
                           .getVertexCoordinates();
                     double[] p3 = facet.getDart().getInvolution(1)
                           .getInvolution(0)
                           .getVertexCoordinates();

                     double[] edge1 = new double[p1.length];
                     double[] edge2 = new double[p1.length];

                     for (int i = 0; i < p1.length; i++)
                     {
                        edge1[i] = p2[i] - p1[i];
                        edge2[i] = p3[i] - p1[i];
                     }

                     double[] facetNormal = Base.cross3(edge1, edge2);

                     double lambda = Base.cutPlaneEdgeV2(facetNormal, p1, pointCoordinates, ray);

                     if (lambda != Double.NaN)
                     {
                        double[] pointOnPlane = new double[p1.length];
                        for (int i = 0; i < p1.length; i++)
                        {
                           pointOnPlane[i] = pointCoordinates[i] + lambda * ray[i];
                        }

                        if (facet.containsPoint(tolerance, pointOnPlane))
                        {
                           nCuts++;
                        }
                     }
                  }

                  if (nCuts % 2 != 0)
                  {
                     return true; // odd number of cuts
                  }
                  else
                     return false; // equal number of cuts

               }
               else
               {
                  System.out.println("Cell.containsPoint: not implemented for dimension = " + dimension + "!");
               }
      return ret;
   }

   /**
    * intersects this cell with another cell and returns resulting cells
    * 
    * @param anotherCell
    * @return List of resulting cells
    */
   public CellList intersect(Cell anotherCell)
   {
      return this.intersect(anotherCell, GMap.TOLERANCE_EPS);
   }

   /**
    * intersects this cell with another cell and returns resulting cells
    * 
    * @param anotherCell
    * @param tolerance
    * @return List of resulting cells
    */
   public CellList intersect(Cell anotherCell, double tolerance)
   {
      CellList retLst = new CellList();

      if (this == anotherCell)
      {
         retLst.setDimension(0);
         retLst.add(this);
      }
      else
         if (this.dimension == 0 && anotherCell.dimension == 0)
         {
            retLst.setDimension(0);

            if (this.containsPoint(tolerance, anotherCell.vertex))
            {
               retLst.add(this);
            }
         }
         else
            if ((this.dimension == 0 && anotherCell.dimension == 1) ||
                  (this.dimension == 1 && anotherCell.dimension == 0))
            {
               // point on edge
            }

      return retLst;
   }

   /**
    * calculates the DualGeometry assuming this cell is an n-Cell
    *
    * @return
    */
   public double[] getDualGeometry()
   {
      return getDualGeometry(true);// or better true
   }

   /**
    * calculates the DualGeometry assuming this cell is an n-Cell
    * 
    * @param exactCalc
    *           for dim R^n and k<=n+1 points use exact calculations
    * @return
    */
   public double[] getDualGeometry(boolean exactCalc)
   {
      // planned to be a replacement for Base.getDualGeometry and Base.getDualGeometryByArr
      List<double[]> points = getVertexCoordinates();

      if (points.size() == 1)
      {
         return points.get(0);
      }
      else
         if (points.size() == 2)
         {
            double[] ret = new double[points.get(0).length];

            double[] p1 = points.get(0);
            double[] p2 = points.get(1);

            for (int i = 0; i < p1.length; i++)
            {
               ret[i] = (p1[i] + p2[i]) * 0.5;
            }

            return ret;
         }
         else
            if (exactCalc && points.size() <= points.get(0).length + 1)
            {
               Sphere sph = Sphere.getCircumscribedSphere(points,true);
               if (sph != null)
                  return sph.getCenter();
               else
                  return getDualGeometry(false);
            }
            else
               if (points.size() == 3)
               {

                  double le1 = Base.distanceSq(points.get(1), points.get(0));
                  double le2 = Base.distanceSq(points.get(2), points.get(1));
                  double le3 = Base.distanceSq(points.get(0), points.get(2));

                  double[] p1 = points.get(0);
                  double[] p2 = points.get(1);
                  double[] p3 = points.get(2);

                  double[] ret = new double[p1.length];

                  if (le1 == le2 && le2 == le3)
                  {
                     // equilateral
                     for (int i = 0; i < p1.length; i++)
                     {
                        ret[i] = (p1[i] + p2[i] + p3[i]) / 3;
                     }
                  }
                  else
                  {
                     boolean isIn = false;
//                     Sphere sph = new Sphere();

                     if (le1 > le2)
                     {
                        if (le1 > le3)
                        {
                           for (int i = 0; i < p1.length; i++)
                           {
                              ret[i] = (p1[i] + p2[i]) * 0.5;
                           }
                           isIn = Base.distanceSq(ret, p3) - (0.5 * le1) < GMap.TOLERANCE_EPS;
                        }
                        else
                        {
                           for (int i = 0; i < p1.length; i++)
                           {
                              ret[i] = (p1[i] + p3[i]) * 0.5;
                           }
                           isIn = Base.distanceSq(ret, p2) - (0.5 * le3) < GMap.TOLERANCE_EPS;
                        }
                     }
                     else
                     {
                        if (le2 > le3)
                        {
                           for (int i = 0; i < p1.length; i++)
                           {
                              ret[i] = (p2[i] + p3[i]) * 0.5;
                           }
                           isIn = Base.distanceSq(ret, p1) - (0.5 * le2) < GMap.TOLERANCE_EPS;
                        }
                        else
                        {
                           for (int i = 0; i < p1.length; i++)
                           {
                              ret[i] = (p1[i] + p3[i]) * 0.5;
                           }
                           isIn = Base.distanceSq(ret, p2) - (0.5 * le3) < GMap.TOLERANCE_EPS;
                        }
                     }

                     if (!isIn)
                     {
                        ret = Sphere.getCircumscribedSphere(getVertexCoordinates(),true).getCenter();
                     }
                  }
                  return ret;
               }
               else
                  if (points.size() == 4)
                  {
                     double[][] pointArr = new double[4][];
                     pointArr[0] = points.get(0);
                     pointArr[1] = points.get(1);
                     pointArr[2] = points.get(2);
                     pointArr[3] = points.get(3);

                     double[] ret = new double[pointArr[0].length];

                     double[][] el = new double[6][3];
                     el[0][0] = Base.distanceSq(points.get(1), points.get(0));
                     el[0][1] = 0;
                     el[0][2] = 1;
                     el[1][0] = Base.distanceSq(points.get(1), points.get(2));
                     el[1][1] = 1;
                     el[1][2] = 2;
                     el[2][0] = Base.distanceSq(points.get(0), points.get(2));
                     el[2][1] = 0;
                     el[2][2] = 2;
                     el[3][0] = Base.distanceSq(points.get(3), points.get(0));
                     el[3][1] = 0;
                     el[3][2] = 3;
                     el[4][0] = Base.distanceSq(points.get(1), points.get(3));
                     el[4][1] = 1;
                     el[4][2] = 3;
                     el[5][0] = Base.distanceSq(points.get(2), points.get(3));
                     el[5][1] = 2;
                     el[5][2] = 3;

                     double max = el[0][0];
                     int id1 = (int) el[0][1];
                     int id2 = (int) el[0][2];

                     for (int i = 1; i < 6; i++)
                     {
                        if (max < el[i][0])
                        {
                           max = el[i][0];
                           id1 = (int) el[i][1];
                           id2 = (int) el[i][2];
                        }
                     }

                     for (int i = 0; i < pointArr[0].length; i++)
                     {
                        ret[i] = (pointArr[id1][i] + pointArr[id2][i]) * 0.5;
                     }

                     boolean isIn = true;
                     int[] otherIds = new int[2];

                     int i = 0;

                     int curOI = 0;
                     while (i < 4)
                     {
                        if (i != id1 && i != id2)
                        {
                           isIn = isIn && Base.distanceSq(ret, pointArr[i]) - (0.5 * max) < GMap.TOLERANCE_EPS;;
                           otherIds[curOI] = i;
                           curOI++;
                        }
                        i++;
                     }

                     if (!isIn)
                     {
                        List<double[]> subal = new ArrayList<double[]>();
                        subal.add(pointArr[id1]);
                        subal.add(pointArr[id2]);
                        subal.add(pointArr[otherIds[0]]);
                        Sphere sph1 = Sphere.getCircumscribedSphere(subal,true);

                        subal.remove(2);
                        subal.add(pointArr[otherIds[1]]);
                        Sphere sph2 = Sphere.getCircumscribedSphere(subal,true);

                        boolean sph1IsValid = sph1 != null && sph1.isIn(pointArr[otherIds[1]]);
                        boolean sph2IsValid = sph2 != null && sph2.isIn(pointArr[otherIds[0]]);

                        if (sph1IsValid && sph2IsValid)
                        {
                           if (sph1.getRadius() < sph2.getRadius())
                           {
                              ret = sph1.getCenter();
                           }
                           else
                           {
                              ret = sph2.getCenter();
                           }
                        }
                        else
                           if (sph1IsValid && !sph2IsValid)
                           {
                              ret = sph1.getCenter();
                           }
                           else
                              if (!sph1IsValid && sph2IsValid)
                              {
                                 ret = sph2.getCenter();
                              }
                              else
                              {
                                 Sphere sph = Sphere.getCircumscribedSphere(points,true);
                                 ret = sph.getCenter();
                              }
                     }

                     return ret;
                  }
                  else
                  {
                     return Sphere.getCircumscribedSphere(getVertexCoordinates(),true).getCenter();
                  }
   }

   public String print()
   {
      String ret = this + "\n";
      ret = ret + "\tdimension : " + dimension + "\n";
      if (dimension == 0)
         ret = ret + "\tcoordinates : " + Arrays.toString(vertex) + "\n";
      if (cellDart != null)
      {
         ret = ret + "\tcell dart:\n";
         ret = ret + cellDart.print();
      }
      else
         ret = ret + "\tcell dart: null";

      if (data != null)
         ret = ret + "\tcell data(" + data.getClass().getName() + "): " + data;

      return ret;
   }

   @Override
   public int getDimension()
   {
      // TODO Auto-generated method stub
      return dimension;
   }

   @Override
   public double getMeasure()
   {
      // TODO Auto-generated method stub
      if (dimension == 0)
         return 0;
      else
         if (dimension == 1)
         {
            // TODO: return edge length
            double[] c1 = cellDart.getIncidentCell(0).vertex;
            double[] c2 = cellDart.getInvolution(0).getIncidentCell(0).vertex;
            // System.out.println(Arrays.toString(c1)+"<>"+Arrays.toString(c2));
            return Base.distance(c1, c2);
         }
         else
            if (dimension == 2)
            {
               // TODO: return facet area

               Dart facetDart = this.cellDart;
               if (facetDart.getIncidentCell(0).getVertex().length == 3)
               {
                  // old version using Base.polygonArea3(vertices);
                  // ArrayList<double[]> facetVertexList = new ArrayList<double[]>();
                  //
                  // Dart d = facetDart;
                  // if (d.getPolarity() == -1)
                  // d = facetDart.getInvolution(1);
                  //
                  // Dart nd = d.getInvolution(0);
                  // while (nd != d.getInvolution(1))
                  // {
                  // facetVertexList.add(nd.getIncidentCell(0).getVertex());
                  // nd = nd.getInvolution(1).getInvolution(0);
                  // }
                  //
                  // double[][] vertices = new double[facetVertexList.size()][];
                  //
                  // facetVertexList.toArray(vertices);
                  // return Base.polygonArea3(vertices);

                  // new version, better do calculations directly here

                  Dart d = facetDart;
                  double[] p1 = d.getVertexCoordinates();
                  double[] p2 = d.getInvolution(0)
                        .getVertexCoordinates();
                  double[] p3 = d.getInvolution(1)
                        .getInvolution(0)
                        .getVertexCoordinates();

                  double[] edge1 = new double[p1.length];
                  double[] edge2 = new double[p1.length];

                  for (int i = 0; i < p1.length; i++)
                  {
                     edge1[i] = p2[i] - p1[i];
                     edge2[i] = p3[i] - p1[i];
                  }

                  double[] facetNormal = Base.cross3(edge1, edge2);
                  double normN = Base.norm(facetNormal);
                  for (int i = 0; i < p1.length; i++)
                  {
                     facetNormal[i] = facetNormal[i] / normN;
                  }

                  double[] area = new double[3];

                  if (d.getPolarity() == -1)
                     d = facetDart.getInvolution(1);

                  Dart nd = d.getInvolution(0);

                  while (nd != d.getInvolution(1))
                  {
                     double[] cp = Base.cross3(nd.getInvolution(0).getVertexCoordinates(), nd.getVertexCoordinates());
                     area[0] = area[0] + cp[0];
                     area[1] = area[1] + cp[1];
                     area[2] = area[2] + cp[2];

                     nd = nd.getInvolution(1).getInvolution(0);
                  }

                  double[] cp = Base.cross3(nd.getInvolution(0).getVertexCoordinates(), nd.getVertexCoordinates());
                  area[0] = area[0] + cp[0];
                  area[1] = area[1] + cp[1];
                  area[2] = area[2] + cp[2];

                  return 0.5 * Math.abs(Base.dot(facetNormal, area));
               }
               else
                  if (facetDart.getIncidentCell(0).getVertex().length == 2)
                  {
                     // old version using Base.polygonArea2(vertices);
                     // ArrayList<double[]> facetVertexList = new ArrayList<double[]>();
                     //
                     // Dart d = facetDart;
                     // if (d.getPolarity() == -1)
                     // d = facetDart.getInvolution(1);
                     //
                     // Dart nd = d.getInvolution(0);
                     // while (nd != d.getInvolution(1))
                     // {
                     // facetVertexList.add(nd.getIncidentCell(0).getVertex());
                     // nd = nd.getInvolution(1).getInvolution(0);
                     // }
                     //
                     // double[][] vertices = new double[facetVertexList.size()][];
                     //
                     // facetVertexList.toArray(vertices);
                     // return Base.polygonArea2(vertices);
                     // new version, better do calculations directly here
                     Dart d = facetDart;

                     double area = 0;

                     if (d.getPolarity() == -1)
                        d = facetDart.getInvolution(1);

                     Dart nd = d.getInvolution(0);
                     while (nd != d.getInvolution(1))
                     {
                        double[] p1 = nd.getInvolution(0).getVertexCoordinates();
                        double[] p2 = nd.getVertexCoordinates();
                        double cp = p1[0] * p2[1] - p1[1] * p2[0];
                        area = area + cp;

                        nd = nd.getInvolution(1).getInvolution(0);
                     }

                     double[] p1 = nd.getInvolution(0).getVertexCoordinates();
                     double[] p2 = nd.getVertexCoordinates();
                     double cp = p1[0] * p2[1] - p1[1] * p2[0];
                     area = area + cp;

                     return 0.5 * Math.abs(area);
                  }
                  else
                     throw new IllegalArgumentException("This operations is implemented only for 2D and 3D points!");

            }
            else
               if (dimension == 3)
               {

                  CellList facets = getBorderingCells(2);

                  double volume = 0;
                  for (int f = 0; f < facets.size(); f++)
                  {
                     Cell facet = facets.get(f);
                     Dart facetDart = facet.getDart();
                     Dart d = facetDart;

                     double[] p1 = d.getVertexCoordinates();
                     double[] p2 = d.getInvolution(0)
                           .getVertexCoordinates();
                     double[] p3 = d.getInvolution(1)
                           .getInvolution(0)
                           .getVertexCoordinates();

                     double[] edge1 = new double[p1.length];
                     double[] edge2 = new double[p1.length];

                     for (int i = 0; i < p1.length; i++)
                     {
                        edge1[i] = p2[i] - p1[i];
                        edge2[i] = p3[i] - p1[i];
                     }

                     double[] facetNormal = Base.cross3(edge1, edge2);
                     double normN = Base.norm(facetNormal);
                     for (int i = 0; i < p1.length; i++)
                     {
                        facetNormal[i] = facetNormal[i] / normN;
                     }

                     double[] area = new double[3];

                     if (d.getPolarity() == -1)
                        d = facetDart.getInvolution(1);

                     Dart nd = d.getInvolution(0);
                     while (nd != d.getInvolution(1))
                     {
                        double[] cp = Base.cross3(nd.getInvolution(0).getVertexCoordinates(),
                              nd.getVertexCoordinates());

                        area[0] = area[0] + cp[0];
                        area[1] = area[1] + cp[1];
                        area[2] = area[2] + cp[2];

                        nd = nd.getInvolution(1).getInvolution(0);
                     }

                     double[] cp = Base.cross3(nd.getInvolution(0).getVertexCoordinates(), nd.getVertexCoordinates());
                     area[0] = area[0] + cp[0];
                     area[1] = area[1] + cp[1];
                     area[2] = area[2] + cp[2];

                     volume = volume + (Base.dot3(facetNormal, p1) * 0.5 * Math.abs(Base.dot(facetNormal, area)));
                  }
                  return (1. / 3.) * Math.abs(volume);
                  // TODO: return polyhedron volume
                  // Base#polyhedronVolume(double[][][])

               }
      return 0;
   }

   @Override
   public CellList getBorderingCells()
   {
      // TODO Auto-generated method stub
      return getBorderingCells(this.dimension - 1);
   }

   @Override
   public CellList getNeighboringCells()
   {
      // TODO Auto-generated method stub
      return getNeighboringCells(this.dimension);
   }

   @Override
   public CellList getBorderingCells(int dim)
   {
      // System.out.println(">>>> getBorderingCells start");
      if (dim < this.dimension && dim >= 0)
      {

         CellList cells = new CellList(dim);
         ArrayList<Dart> iOrbit = cellDart.getIOrbit(dimension);

         while (!iOrbit.isEmpty())
         {
            ArrayList<Dart> vDarts = iOrbit.get(0).getIOrbit(dim);

            cells.add(vDarts.get(0).getIncidentCell(dim));
            iOrbit.removeAll(vDarts);
         }

         // System.out.println(">>>>"+cells.getList());
         return cells;
      }
      return null;
   }

   @Override
   public CellList getNeighboringCells(int dim)
   {
      // TODO Auto-generated method stub
      Set<Cell> neighborCells = new HashSet<Cell>();

      ArrayList<Dart> iOrbit = cellDart.getIOrbit(dimension);

      for (int i = 0; i < iOrbit.size(); i++)
      {
         neighborCells.add(iOrbit.get(i).getInvolution(dimension).getIncidentCell(dim));
      }

      return new CellList(neighborCells);
   }

   @Override
   public CellList getVertices()
   {
      // TODO: may be not efficient enough ...
      CellList vCells = new CellList(0);
      ArrayList<Dart> iOrbit = cellDart.getIOrbit(dimension);

      while (!iOrbit.isEmpty())
      {
         ArrayList<Dart> vDarts = iOrbit.get(0).getIOrbit(0);

         vCells.add(vDarts.get(0).getIncidentCell(0));
         iOrbit.removeAll(vDarts);
      }

      return vCells;
   }

   @Override
   public List<double[]> getVertexCoordinates()
   {
      List<double[]> vrtxCoords = new ArrayList<double[]>();
      ArrayList<Dart> iOrbit = cellDart.getIOrbit(dimension);

      while (!iOrbit.isEmpty())
      {
         ArrayList<Dart> vDarts = iOrbit.get(0).getIOrbit(0);

         vrtxCoords.add(vDarts.get(0).getVertexCoordinates());
         iOrbit.removeAll(vDarts);
      }

      return vrtxCoords;
   }

   @Override
   public Sphere getCircumscribingCircle() throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public double[] getCenterOfMass()
   {
      if (dimension == 0)
      {
         return vertex;
      }
      else
         if (dimension == 1)
         {
            double[] v1 = cellDart.getVertexCoordinates();
            double[] v2 = cellDart.getInvolution(0).getVertexCoordinates();

            double[] ret = new double[v1.length];
            for (int i = 0; i < v1.length; i++)
            {
               ret[i] = (v1[i] + v2[i]) * 0.5;
            }
            return ret;
         }
         else
         {
            return getVertices().getCenterOfMass();
         }
   }

}
