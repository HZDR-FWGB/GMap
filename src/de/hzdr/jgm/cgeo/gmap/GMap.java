package de.hzdr.jgm.cgeo.gmap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.IdentityHashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import de.hzdr.jgm.cgeo.Base;
import de.hzdr.jgm.cgeo.NearestNeighbor;

/**
 * Class for Objects, represented by a GMap
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class GMap implements GGeometry
{

   // static utility
   /**
    * tolerance
    */
   static final double TOLERANCE_EPS = 0.0000000001;

   static public double getInternalTolerance()
   {
      return TOLERANCE_EPS;
   }
   // static creation

   /**
    * Creates a single quad based on two points (axis aligned) or a sequence of four points.
    * 
    * @param points
    * @return
    */
   public static GMap createQuad(double[]... points)
   {
      if (points.length == 2)
      {
         double[][] pointSequence;

         if (points[0].length == 2)
         {
            pointSequence = new double[4][2];
            pointSequence[0] = points[0];
            pointSequence[1][0] = points[1][0];
            pointSequence[1][1] = points[0][1];
            pointSequence[2] = points[1];
            pointSequence[3][0] = points[0][0];
            pointSequence[3][1] = points[1][1];
         }
         else
            if (points[0].length == 3)
            {
               pointSequence = new double[4][3];

               pointSequence[0] = points[0];
               pointSequence[1][0] = points[1][0];
               pointSequence[1][1] = points[0][1];
               pointSequence[1][2] = points[0][2];
               pointSequence[2] = points[1];
               pointSequence[3][0] = points[0][0];
               pointSequence[3][1] = points[1][1];
               pointSequence[3][2] = points[1][2];
            }
            else
            {
               throw new IllegalArgumentException("operation only allowed for dimensions 2 and 3");
            }

         for (int i = 0; i < 4; i++)
            System.out.println(">>" + Arrays.toString(pointSequence[i]));
         return createPolygon(pointSequence);
      }
      else
         return createPolygon(points);
   }

   /**
    * Creates a single polygon based on a sequence of points.
    * 
    * @param pointSequence
    * @return
    */
   public static GMap createPolygon(double[]... pointSequence)
   {
      GMap gmo = new GMap(2);

      CellList vrtxCells = gmo.addPoints(pointSequence);
      CellList edgeCells = new CellList(1);
      for (int e = 0; e < pointSequence.length; e++)
      {
         Cell v1 = vrtxCells.get(e);
         Cell v2;
         if (e < pointSequence.length - 1)
            v2 = vrtxCells.get(e + 1);
         else
            v2 = vrtxCells.get(0);

         edgeCells.add(gmo.addCell(v1, v2));;
      }
      gmo.addCell(edgeCells);
      gmo.connectCells();
      return gmo;
   }

   /**
    * creates a 2-GMap representing a polygonal mesh
    * 
    * Note use this instead of {@link #loadVTK(String) loadVTK(String)} for PolyData with triangles/polygons
    * 
    * @param pointList
    *           points
    * @param polygonList
    *           polygon definitions
    * @return new GMap object
    */
   @SuppressWarnings("unchecked")
   public static GMap createPolygonMesh(List<double[]> pointList, List<int[]> polygonList)
   {
      // step 1: building incidence graph
      long start = System.currentTimeMillis();
      int dimension = 2;
      Map<Cell, Set<Cell>[]> incidenceGraph = new IdentityHashMap<Cell, Set<Cell>[]>();

      GMap gmo = new GMap(dimension);

      if (pointList.size() > 0 && polygonList.size() > 0)
      {
         // building edge-vertex list and face-edge list

         // create all vertex cells
         ArrayList<Cell> vertices = new ArrayList<Cell>();

         ArrayList<Integer>[] edgesAtVertex = (ArrayList<Integer>[]) new ArrayList[pointList.size()];
         ArrayList<Integer>[] facetsAtVertex = (ArrayList<Integer>[]) new ArrayList[pointList.size()];

         for (int p = 0; p < pointList.size(); p++)
         {
            Cell v = new Cell(pointList.get(p));
            vertices.add(v);
            Set<Cell>[] incidentCells = (HashSet<Cell>[]) new HashSet[2];
            incidentCells[0] = new HashSet<Cell>();
            incidentCells[1] = new HashSet<Cell>();
            incidenceGraph.put(v, incidentCells);

            edgesAtVertex[p] = new ArrayList<Integer>();
            facetsAtVertex[p] = new ArrayList<Integer>();
         }

         ArrayList<Cell> edges = new ArrayList<Cell>();
         ArrayList<Cell> facets = new ArrayList<Cell>();
         for (int p = 0; p < polygonList.size(); p++)
         {
            Cell newFacet = new Cell(2);
            facets.add(newFacet);
            Set<Cell>[] incidentCells = (HashSet<Cell>[]) new HashSet[2];
            incidentCells[0] = new HashSet<Cell>();
            incidentCells[1] = new HashSet<Cell>();
            incidenceGraph.put(newFacet, incidentCells);

            int[] polygonIds = polygonList.get(p);

            for (int e = 0; e < polygonIds.length; e++)
            {
               int vid1 = polygonIds[e];
               facetsAtVertex[vid1].add(p);
               int vid2;
               if (e < polygonIds.length - 1)
                  vid2 = polygonIds[e + 1];
               else
                  vid2 = polygonIds[0];
               facetsAtVertex[vid2].add(p);

               ArrayList<Integer> incidentEdgeIds1 = edgesAtVertex[vid1];
               ArrayList<Integer> incidentEdgeIds2 = edgesAtVertex[vid2];

               int tmpEdgeID = -1;

               for (int e1 = 0; e1 < incidentEdgeIds1.size(); e1 = e1 + 1)
               {
                  Integer egdeId = incidentEdgeIds1.get(e1);
                  for (int e11 = 0; e11 < incidentEdgeIds2.size(); e11 = e11 + 1)
                  {
                     Integer egdeId2 = incidentEdgeIds2.get(e11);
                     if (egdeId == egdeId2)
                     {
                        // edge between both vertices already existing
                        tmpEdgeID = egdeId;
                        break;
                     }
                  }
                  if (tmpEdgeID != -1)
                     break;
               }

               Cell edge;
               if (tmpEdgeID != -1)
               {
                  edge = edges.get(tmpEdgeID);
                  incidenceGraph.get(edge)[1].add(newFacet);
               }
               else
               {
                  edge = new Cell(1);

                  Integer nEID = edges.size();

                  edgesAtVertex[vid1].add(nEID);
                  edgesAtVertex[vid2].add(nEID);

                  edges.add(edge);

                  Set<Cell>[] incidentCellsE = (HashSet<Cell>[]) new HashSet[2];
                  incidentCellsE[0] = new HashSet<Cell>();
                  incidentCellsE[1] = new HashSet<Cell>();

                  incidentCellsE[0].add(vertices.get(vid1));
                  incidentCellsE[0].add(vertices.get(vid2));
                  incidentCellsE[1].add(newFacet);

                  incidenceGraph.put(edge, incidentCellsE);

                  incidenceGraph.get(vertices.get(vid1))[1].add(edge);
                  incidenceGraph.get(vertices.get(vid2))[1].add(edge);
               }

               incidentCells[0].add(edge);
            }
         }
         System.out.println("Building incidence graph for polygonal mesh ( "
               + (System.currentTimeMillis() - start) + "ms)");

         gmo.cells[0].addAll(vertices);
         gmo.cells[1].addAll(edges);
         gmo.cells[2].addAll(facets);

         Map<Cell, ArrayList<Dart>> cellDarts = new IdentityHashMap<Cell, ArrayList<Dart>>();
         ArrayList<Dart> newDartLst = new ArrayList<Dart>();

         for (int c = 0; c < facets.size(); c += 1)
         {
            newDartLst.addAll(gmo.addNCell(facets.get(c), incidenceGraph, cellDarts));
         }
         // gmo.check();
      }

      System.out.println("Building GMap for polygonal mesh ( "
            + (System.currentTimeMillis() - start) + "ms)");

      return gmo;
   }

   /**
    * creates a 3-GMap representing a tetrahedral mesh
    * 
    * Notes:
    * 
    * 1)tetrahedron definitions according to VTK convention
    * 
    * 2) use this instead of {@link #loadVTK(String) loadVTK(String)} for unstructured grids with cell type 10
    * 
    * @param pointList
    *           points
    * @param tetrahedronList
    *           polygon definitions
    * @return new GMap object
    */
   @SuppressWarnings("unchecked")
   public static GMap createTetrahedronMesh(List<double[]> pointList, List<int[]> tetrahedronList)
   {
      // step 1: building incidence graph
      long start = System.currentTimeMillis();
      int dimension = 3;
      Map<Cell, Set<Cell>[]> incidenceGraph = new IdentityHashMap<Cell, Set<Cell>[]>();

      GMap gmo = new GMap(dimension);
      if (pointList.size() > 0 && tetrahedronList.size() > 0)
      {
         // building edge-vertex list and face-edge list

         // create all vertex cells
         ArrayList<Cell> vertices = new ArrayList<Cell>();

         ArrayList<Integer>[] edgesAtVertex = (ArrayList<Integer>[]) new ArrayList[pointList.size()];
         ArrayList<Integer>[] facetsAtVertex = (ArrayList<Integer>[]) new ArrayList[pointList.size()];
         ArrayList<Integer>[] cellsAtVertex = (ArrayList<Integer>[]) new ArrayList[pointList.size()];

         for (int p = 0; p < pointList.size(); p++)
         {
            Cell v = new Cell(pointList.get(p));
            vertices.add(v);
            Set<Cell>[] incidentCells = (HashSet<Cell>[]) new HashSet[2];
            incidentCells[0] = new HashSet<Cell>();
            incidentCells[1] = new HashSet<Cell>();
            incidenceGraph.put(v, incidentCells);

            edgesAtVertex[p] = new ArrayList<Integer>();
            facetsAtVertex[p] = new ArrayList<Integer>();
            cellsAtVertex[p] = new ArrayList<Integer>();
         }

         ArrayList<Cell> edges = new ArrayList<Cell>();
         ArrayList<Cell> facets = new ArrayList<Cell>();
         ArrayList<Cell> tetras = new ArrayList<Cell>();

         for (int p = 0; p < tetrahedronList.size(); p++)
         {
            Cell newTetra = new Cell(3);
            tetras.add(newTetra);
            Set<Cell>[] incidentCells = (HashSet<Cell>[]) new HashSet[2];
            incidentCells[0] = new HashSet<Cell>();
            incidentCells[1] = new HashSet<Cell>();
            incidenceGraph.put(newTetra, incidentCells);

            int[] tetraIds = tetrahedronList.get(p);

            for (int i = 0; i < tetraIds.length; i++)
               cellsAtVertex[tetraIds[i]].add(p);

            // corresponds to VTK cell type 10 (tetrahedron)
            // here the cell type is defined - adapt for other 3D meshes with equal cell types
            int[][] tetraFaceIds =
            {
                  { tetraIds[0], tetraIds[1], tetraIds[2] },
                  { tetraIds[3], tetraIds[1], tetraIds[0] },
                  { tetraIds[3], tetraIds[2], tetraIds[1] },
                  { tetraIds[3], tetraIds[0], tetraIds[2] } };

            for (int cf = 0; cf < tetraFaceIds.length; cf += 1)
            {
               int[] face = tetraFaceIds[cf];

               // System.out.println("Facet: "+Arrays.toString(face));
               Cell newFacet;
               // new or existing
               int tmpFaceId = -1;
               boolean isNew = false;
               ArrayList<Integer>[] incidentFacets = (ArrayList<Integer>[]) new ArrayList[face.length];
               for (int i = 0; i < face.length; i += 1)
               {
                  incidentFacets[i] = facetsAtVertex[face[i]];
                  // System.out.println("incident facets: "+incidentFacets[i]);
                  if (incidentFacets[i].size() == 0)
                  {
                     // System.out.println("... new because of unused vertex");
                     isNew = true;
                     break;
                  }
               }

               if (!isNew)
               {
                  // find one id in all lists
                  for (int i = 0; i < incidentFacets[0].size(); i++)
                  {
                     Integer id = incidentFacets[0].get(i);
                     int j = 1;
                     while (j < incidentFacets.length && incidentFacets[j].contains(id))
                        j++;

                     if (j == incidentFacets.length)
                     {
                        tmpFaceId = id;
                        break;
                     }
                  }

               }

               Set<Cell>[] incidentCellsF = (HashSet<Cell>[]) new HashSet[2];
               if (tmpFaceId > -1)
               {
                  // System.out.println("... already existing");
                  newFacet = facets.get(tmpFaceId);

                  incidentCellsF = incidenceGraph.get(newFacet);
                  incidentCellsF[1].add(newTetra);
               }
               else
               {
                  // System.out.println("... new");
                  newFacet = new Cell(2);
                  Integer nFID = facets.size();

                  for (int i = 0; i < face.length; i += 1)
                  {
                     // System.out.println("... at vrtx "+face[i]+" add facet "+nFID);
                     facetsAtVertex[face[i]].add(nFID);
                     // System.out.println("...."+facetsAtVertex[face[i]]);
                  }

                  facets.add(newFacet);

                  incidentCellsF[0] = new HashSet<Cell>();
                  incidentCellsF[1] = new HashSet<Cell>();
                  incidentCellsF[1].add(newTetra);

                  incidenceGraph.put(newFacet, incidentCellsF);
               }

               for (int e = 0; e < face.length; e++)
               {
                  int vid1 = face[e];
                  // facetsAtVertex[vid1].add(p);
                  int vid2;
                  if (e < face.length - 1)
                     vid2 = face[e + 1];
                  else
                     vid2 = face[0];
                  // facetsAtVertex[vid2].add(p);

                  ArrayList<Integer> incidentEdgeIds1 = edgesAtVertex[vid1];
                  ArrayList<Integer> incidentEdgeIds2 = edgesAtVertex[vid2];

                  int tmpEdgeID = -1;

                  for (int e1 = 0; e1 < incidentEdgeIds1.size(); e1 = e1 + 1)
                  {
                     Integer egdeId = incidentEdgeIds1.get(e1);
                     for (int e11 = 0; e11 < incidentEdgeIds2.size(); e11 = e11 + 1)
                     {
                        Integer egdeId2 = incidentEdgeIds2.get(e11);
                        if (egdeId == egdeId2)
                        {
                           // edge between both vertices already existing
                           tmpEdgeID = egdeId;
                           break;
                        }
                     }
                     if (tmpEdgeID != -1)
                        break;
                  }

                  Cell edge;
                  if (tmpEdgeID != -1)
                  {
                     edge = edges.get(tmpEdgeID);
                     incidenceGraph.get(edge)[1].add(newFacet);
                  }
                  else
                  {
                     edge = new Cell(1);

                     Integer nEID = edges.size();

                     edgesAtVertex[vid1].add(nEID);
                     edgesAtVertex[vid2].add(nEID);

                     edges.add(edge);

                     Set<Cell>[] incidentCellsE = (HashSet<Cell>[]) new HashSet[2];
                     incidentCellsE[0] = new HashSet<Cell>();
                     incidentCellsE[1] = new HashSet<Cell>();

                     incidentCellsE[0].add(vertices.get(vid1));
                     incidentCellsE[0].add(vertices.get(vid2));
                     incidentCellsE[1].add(newFacet);

                     incidenceGraph.put(edge, incidentCellsE);

                     incidenceGraph.get(vertices.get(vid1))[1].add(edge);
                     incidenceGraph.get(vertices.get(vid2))[1].add(edge);
                  }

                  incidentCellsF[0].add(edge);
               }
               incidentCells[0].add(newFacet);
            }
         }

         System.out.println("Building incidence graph for tetrahedral mesh ( "
               + (System.currentTimeMillis() - start) + "ms)");

         gmo.cells[0].addAll(vertices);
         gmo.cells[1].addAll(edges);
         gmo.cells[2].addAll(facets);
         gmo.cells[3].addAll(tetras);

         // System.out.println("Vertices: "+gmo.cells[0].size()+", edge: "+gmo.cells[1].size()+"; facets:
         // "+gmo.cells[2].size()+"; polyhedrons: "+gmo.cells[3].size());

         Map<Cell, ArrayList<Dart>> cellDarts = new IdentityHashMap<Cell, ArrayList<Dart>>();
         ArrayList<Dart> newDartLst = new ArrayList<Dart>();

         for (int c = 0; c < tetras.size(); c += 1)
         {
            newDartLst.addAll(gmo.addNCell(tetras.get(c), incidenceGraph, cellDarts));
         }
      }
      System.out.println("Building GMap for tetrahedral mesh ( " + (System.currentTimeMillis() - start) + "ms)");
      return gmo;

   }

   /**
    * creats a new GMap representing the Delaunay tessellation of the given coordinates
    * 
    * @param points
    *           set of coordinates
    * @return
    * @throws NotImplementedException
    */
   public static GMap createDelaunayTessellation(double[]... points) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   /**
    * creats a new GMap representing the Voronoi tessellation of the given coordinates
    * 
    * @param points
    *           set of coordinates
    * @return
    * @throws NotImplementedException
    */
   public static GMap createVoronoiTessellation(double[]... points) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   /**
    * creats a new GMap representing the convex hull (simplicial complex) of the given coordinates
    * 
    * @param points
    *           set of coordinates
    * @return
    * @throws NotImplementedException
    */
   public static GMap createConvexHull(double[]... points) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   /**
    * creats a new GMap representing the lower convex hull (simplicial complex) of the given coordinates
    * 
    * @param points
    *           set of coordinates
    * @return
    * @throws NotImplementedException
    */
   public static GMap createLowerCovexHull(double[]... points) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   static Map<Cell, Cell> createConnectivityFunction(GMap gmap1, GMap gmap2)
   {
      CellList mapVertices1 = gmap1.getVertexCells();
      CellList mapVertices2 = gmap2.getVertexCells();

      // connectivity function
      // association between on vertex v_i to the n-cell c_j incident to its nearest neighbor v_i is in c_j
      Map<Cell, Cell> connectivityFunction = new IdentityHashMap<Cell, Cell>();

      double[][] vertexArr1 = new double[mapVertices1.size()][];
      for (int i = 0; i < mapVertices1.size(); i++)
      {
         vertexArr1[i] = mapVertices1.get(i).getVertex();
      }

      double[][] vertexArr2 = new double[mapVertices2.size()][];
      for (int i = 0; i < mapVertices2.size(); i++)
      {
         vertexArr2[i] = mapVertices2.get(i).getVertex();
      }

      NearestNeighbor nn_Map1 = new NearestNeighbor(vertexArr1);
      NearestNeighbor nn_Map2 = new NearestNeighbor(vertexArr2);

      if ((gmap1.maxDimension == 2 || gmap1.maxDimension == 3) && (gmap2.maxDimension == 2 || gmap2.maxDimension == 3))
      {
         for (int i = 0; i < mapVertices1.size(); i++)
         {
            if (gmap2.cells[gmap2.maxDimension].size() > 1)
            {
               double[] vrtxCoord = vertexArr1[i];
               int id = nn_Map2.getIndex(vrtxCoord);

               CellList incidentNCells = mapVertices2.get(id).getNeighboringCells(gmap2.maxDimension);

               int j = 0;
               while (j < incidentNCells.size() && !incidentNCells.get(j).containsPoint(vrtxCoord))
               {
                  j++;
               }

               if (j < incidentNCells.size())
               {
                  connectivityFunction.put(mapVertices1.get(i), incidentNCells.get(j));
               }
            }
            else
            {
               Cell ncell = gmap2.cells[gmap2.maxDimension].get(0);
               for (int j = 0; j < mapVertices1.size(); j++)
               {
                  double[] cellVrtxCoord = vertexArr1[i];
                  Cell vc = mapVertices2.containsPoint(vertexArr1[i]);
                  if (ncell.containsPoint(cellVrtxCoord))
                  {
                     if (vc != null)
                     {
                        connectivityFunction.put(vc, ncell);
                     }

                  }
                  else
                  {
                     int id = nn_Map2.getIndex(cellVrtxCoord);
                     connectivityFunction.put(vc, mapVertices2.get(id));
                  }
               }
            }

            for (int k = 0; k < mapVertices2.size(); k++)
            {
               if (gmap1.cells[gmap1.maxDimension].size() > 1)
               {
                  double[] vrtxCoord = vertexArr2[k];
                  int id = nn_Map1.getIndex(vrtxCoord);

                  CellList incidentNCells = mapVertices1.get(id).getNeighboringCells(gmap2.maxDimension);

                  int j = 0;
                  while (j < incidentNCells.size() && !incidentNCells.get(j).containsPoint(vrtxCoord))
                  {
                     j++;
                  }

                  if (j < incidentNCells.size())
                  {
                     connectivityFunction.put(mapVertices2.get(k), incidentNCells.get(j));
                  }
               }
               else
               {
                  Cell ncell = gmap1.cells[gmap1.maxDimension].get(0);
                  for (int j = 0; j < mapVertices2.size(); j++)
                  {
                     double[] cellVrtxCoord = vertexArr2[j];
                     Cell vc = mapVertices1.containsPoint(cellVrtxCoord);
                     if (ncell.containsPoint(cellVrtxCoord))
                     {
                        if (vc != null)
                        {
                           connectivityFunction.put(vc, ncell);
                        }

                     }
                     else
                     {
                        int id = nn_Map1.getIndex(cellVrtxCoord);
                        connectivityFunction.put(vc, mapVertices1.get(id));
                     }
                  }
               }
            }
         }
      }
      else
         throw new IllegalArgumentException("This operation is currently implemented for 2D and 3D points only.");

      return connectivityFunction;
   }
   // non static

   private int                    maxDimension            = 3;
   private ArrayList<Dart>        darts                   = new ArrayList<Dart>();
   private ArrayList<Cell>[]      cells;

   private Map<Cell, Set<Cell>[]> temporaryIncidenceGraph = null;

   public GMap()
   {
      // TODO Auto-generated constructor stub
      initCellList();
   }

   public GMap(int maxDim)
   {
      maxDimension = maxDim;
      initCellList();
   }

   // low level interface
   /**
    * gets the list of all darts
    * 
    * @return
    */
   public ArrayList<Dart> getDartList()
   {
      return darts;
   }

   /**
    * resets the complete dart list and build a new cell list
    * 
    * @param darts
    */
   public void setDartList(ArrayList<Dart> newDarts)
   {
      this.darts = new ArrayList<Dart>(darts);
      initCellList();

      for (int d = 0; d < darts.size(); d++)
      {
         Dart dart = this.darts.get(d);
         for (int i = 0; i <= dart.getDimension(); i++)
         {
            if (!cells[i].contains(dart.getIncidentCell(i)))
            {
               if (dart.getIncidentCell(i).getDart() == null)
                  dart.getIncidentCell(i).setDart(dart);
               cells[i].add(dart.getIncidentCell(i));
            }
         }
      }
   }

   /**
    * @return the maxDimension
    */
   public int getMaxDimension()
   {
      return maxDimension;
   }

   /**
    * @param maxDimension
    *           the maxDimension to set
    */
   public void setMaxDimension(int maxDimension)
   {
      this.maxDimension = maxDimension;
   }

   /**
    * @return the temporaryIncidenceGraph
    */
   protected Map<Cell, Set<Cell>[]> getTemporaryIncidenceGraph()
   {
      return temporaryIncidenceGraph;
   }

   /**
    * @param temporaryIncidenceGraph
    *           the temporaryIncidenceGraph to set
    */
   protected void setTemporaryIncidenceGraph(Map<Cell, Set<Cell>[]> temporaryIncidenceGraph)
   {
      this.temporaryIncidenceGraph = temporaryIncidenceGraph;
   }

   /**
    * @return the darts
    */
   protected ArrayList<Dart> getDarts()
   {
      return darts;
   }

   /**
    * @param darts
    *           the darts to set
    */
   protected void setDarts(ArrayList<Dart> darts)
   {
      this.darts = darts;
   }

   /**
    * @return the cells
    */
   public ArrayList<Cell>[] getCells()
   {
      return cells;
   }

   /**
    * @param cells
    *           the cells to set
    */
   protected void setCells(ArrayList<Cell>[] cells)
   {
      this.cells = cells;
   }

   /**
    * adds a new dart to the object without sewing
    * 
    * @param dart
    *           new dart
    * 
    */
   public void addDart(Dart dart)
   {
      addDart(dart, false);
   }

   /**
    * adds a new dart to the object, with optional automated sewing
    * 
    * Note: searching is not fully optimized
    * 
    * @param dart
    *           new dart
    * @param doSew
    *           sew this dart to other darts (inefficient)
    * 
    */
   public void addDart(Dart dart, boolean doSew)
   {

      if (doSew && dart.isOrphan() > -1)
      {
         boolean sewedIn[] = new boolean[dart.getDimension() + 1];
         for (int i = 0; i <= dart.getDimension(); i++)
            sewedIn[i] = dart != dart.getInvolution(i);
         // System.out.println("Add "+ dart.print());
         for (int d = 0; d < darts.size(); d++)
         {
            Dart tD = darts.get(d);
            if (tD.isOrphan() > -1)
            {
               // System.out.println(">to "+ tD.print()+"?");

               boolean isSewed = true;
               for (int i = 0; i <= dart.getDimension(); i++)
               {
                  if (!sewedIn[i] && tD == tD.getInvolution(i) && dart.isInvolution(tD, i))
                  {
                     // System.out.println(">>sew in "+i);
                     dart.sew(tD, i);
                     sewedIn[i] = true;
                     break;
                  }
                  isSewed = isSewed && sewedIn[i];
               }
               if (isSewed)
                  break;
            }
         }
      }

      darts.add(dart);
      for (int i = 0; i <= dart.getDimension(); i++)
      {
         if (!cells[i].contains(dart.getIncidentCell(i)))
         {
            dart.getIncidentCell(i).setDart(dart);
            cells[i].add(dart.getIncidentCell(i));
         }
      }
   }

   /**
    * computes the dual GMap
    * 
    * @return new GMap representing the dual topology/geometry
    */
   public GMap computeDualMap()
   {
      Map<Dart, Dart> dualDartMap = new IdentityHashMap<Dart, Dart>();
      Map<Cell, Cell> dualCellMap = new IdentityHashMap<Cell, Cell>();

      ArrayList<Dart> dualDarts = new ArrayList<Dart>();

      @SuppressWarnings("unchecked")
      ArrayList<Cell>[] dualCells = (ArrayList<Cell>[]) new ArrayList[maxDimension + 1];
      for (int i = 0; i <= maxDimension; i++)
      {
         dualCells[i] = new ArrayList<Cell>();
      }

      // System.out.println("hallo?"+maxDimension);

      // TODO: create dual copies of all embeddings
      for (int i = 0; i <= maxDimension; i++)
      {
         ArrayList<Cell> cList = cells[i];
         ArrayList<Cell> dCList = dualCells[maxDimension - i];
         for (int j = 0; j < cList.size(); j++)
         {
            Cell cell = cList.get(j);
            Cell dualCell = new Cell(maxDimension - cell.getDimension());

            if (cell.getDimension() == maxDimension)
            {
               // following is not implemented yet
               // dualCell.setVertex(Sphere.getCircumscribedSphere(cell.getVertexCoordinates()).getCenter());
               dualCell.setVertex(cell.getDualGeometry());
               // old solution, should not be used
               // dualCell.setVertex(Base.getDualGeometryByArr(
               // cell.getVertexCoordinates(),
               // false, true).getCenterAsArray());
            }

            dualCellMap.put(cell, dualCell);
            dCList.add(dualCell);
         }
      }

      // TODO: create dual copies of all darts, without involutions
      for (int i = 0; i < darts.size(); i++)
      {
         Dart dart = darts.get(i);
         Dart dualDart = new Dart(dart.getDimension());

         for (int j = 0; j <= dart.getDimension(); j++)
         {
            Cell dualCell = dualCellMap.get(dart.getIncidentCell(j));
            // System.out.println(dart.getIncidentCell(j) + "<>"+dualCell);
            if (dualCell.getDart() == null)
               dualCell.setDart(dualDart);

            dualDart.setIncidentCell(dualCell);
         }

         dualDartMap.put(dart, dualDart);
         dualDarts.add(dualDart);
      }

      // TODO: infer involutions
      for (int i = 0; i < darts.size(); i++)
      {
         Dart dart = darts.get(i);
         Dart dualDart = dualDartMap.get(dart);

         for (int j = 0; j <= dart.getDimension(); j++)
         {
            dualDart.getInvolutions()[dualDart.getDimension() - j] = dualDartMap.get(dart.getInvolution(j));
         }
      }

      GMap dualMap = new GMap(maxDimension);
      dualMap.darts = dualDarts;
      dualMap.cells = dualCells;
      return dualMap;
   }

   public boolean check()
   {
      return check(this.darts);
   }

   public boolean check(ArrayList<Dart> drtLst)
   {
      boolean correct = true;
      boolean ret = true;
      Iterator<Dart> dIt = drtLst.iterator();
      int nNonBorderDarts = 0;
      int nBorderdarts = 0;
      int nLowOrphans = 0;

      while (dIt.hasNext())
      {
         Dart d = dIt.next();
         int orphanDim = d.isOrphan();
         if (orphanDim > -1)
         {
            if (orphanDim < this.maxDimension)
            {
               correct = false;
               nLowOrphans++;
            }
            else
            {
               nBorderdarts++;
            }
         }
         else
            nNonBorderDarts++;
      }

      if (!correct)
      {
         System.out.println("GMap.check() Failed! Check for low orphan darts (" + nLowOrphans
               + " detected) - (normal darts: " + nNonBorderDarts + ", bordering darts: " + nBorderdarts + ") -> "
               + (nLowOrphans + nBorderdarts + nNonBorderDarts));
         // return correct;
      }
      else
         System.out.println("GMap.check(): Check for orphan darts (" + nBorderdarts
               + " detected) - (normal darts: " + nNonBorderDarts + ") -> " + (nBorderdarts + nNonBorderDarts)
               + ", passed.");

      if (!correct)
      {
         ret = correct;
         correct = true;
      }

      dIt = drtLst.iterator();
      int fixedPoints1 = 0;
      int fixedPoints2 = 0;
      while (dIt.hasNext())
      {
         Dart d = dIt.next();

         for (int i = 0; i <= maxDimension; i++)
         {
            for (int k = 0; k <= maxDimension; k++)
            {
               int j = i + 2 + k;
               if (j <= maxDimension)
               {
                  if (d == d.getInvolution(i).getInvolution(j))
                  {
                     correct = false;
                     fixedPoints1++;
                  }

                  if (d != d.getInvolution(i).getInvolution(j).getInvolution(i).getInvolution(j))
                  {
                     correct = false;
                     fixedPoints2++;
                  }
               }
            }
         }
      }

      if (!correct)
      {
         System.out.println("GMap.check() Failed check for involution constraints:\n"
               + "d.involution(i).involution(i+2+k)==d (occurance " + fixedPoints1 + " times)\n"
               + "d.involution(i).involution(i+2+k).involution(i).involution(i+2+k)!=d (occurance " + fixedPoints2
               + " times)!");
      }
      else
      {
         System.out.println("GMap.check(): Check for involution constraints passed.");
         if (!ret)
            System.out.println("\t... may be dual GMap ...");
      }

      if (!correct)
      {
         ret = correct;
         correct = true;
      }

      return ret;

   }

   // high level interface implementations
   @Override
   public Cell searchPoint(double... pointCoordinates)
   {
      // TODO Auto-generated method stub
      return searchPoint(GMap.TOLERANCE_EPS, pointCoordinates);
   }

   @Override
   public Cell searchPoint(double tolerance, double... pointCoordinates)
   {

      Cell vrtxC = null;
      Iterator<Cell> iterator = this.cells[0].iterator();
      while (iterator.hasNext() && vrtxC == null)
      {
         Cell emb = iterator.next();
         double[] arr = emb.getVertex();

         if (Base.isEqualCoords(pointCoordinates, arr, tolerance))
         {
            vrtxC = emb;
         }
      }

      return vrtxC;
   }

   @Override
   public Cell getPoint(double... pointCoordinates) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public Cell getPoint(double tolerance, double... pointCoordinates) throws NotImplementedException
   {
      // TODO Auto-generated method stub
      throw new NotImplementedException("Not implemented yet!");
   }

   @SuppressWarnings("unchecked")
   @Override
   public Cell addPoint(double... pointCoordinates)
   {
      Cell vrtxC = new Cell(pointCoordinates);
      cells[0].add(vrtxC);

      // incidence graph block===
      ensureTIG();
      Set<Cell>[] incidentCells = (HashSet<Cell>[]) new HashSet[2];
      incidentCells[0] = new HashSet<Cell>();
      incidentCells[1] = new HashSet<Cell>();
      temporaryIncidenceGraph.put(vrtxC, incidentCells);
      // incidence graph block===

      return vrtxC;
   }

   @Override
   public CellList addPoints(double[]... pointCoordinates)
   {
      CellList vrtxsCell = new CellList(0);
      for (int i = 0; i < pointCoordinates.length; i++)
         vrtxsCell.add(this.addPoint(pointCoordinates[i]));
      return vrtxsCell;
   }

   @Override
   public CellList addPoints(List<double[]> pointCoordinates)
   {
      CellList vrtxsCell = new CellList(0);
      for (int i = 0; i < pointCoordinates.size(); i++)
         vrtxsCell.add(this.addPoint(pointCoordinates.get(i)));
      return vrtxsCell;
   }

   @Override
   public void addVertexList(List<Cell> vertexCells) throws NotImplementedException
   {
      // TODO Auto-generated method stub
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public Cell searchCell(double... pointCoordinates) throws NotImplementedException
   {
      // TODO Auto-generated method stub
      return searchCell(TOLERANCE_EPS, pointCoordinates);
   }

   @Override
   public Cell searchCell(double tolerance, double... pointCoordinates) throws NotImplementedException
   {
      // TODO Auto-generated method stub
      throw new NotImplementedException("Not implemented yet!");
   }

   @SuppressWarnings("unchecked")
   @Override
   public Cell addCell(Cell... borderingCells)
   {
      // TODO Auto-generated method stub
      Cell cell = new Cell(borderingCells[0].getDimension() + 1);
      cell.setDart(borderingCells[0].getDart());

      cells[cell.getDimension()].add(cell);

      // incidence graph block===
      ensureTIG();
      Set<Cell>[] incidentCells = (HashSet<Cell>[]) new HashSet[2];
      incidentCells[0] = new HashSet<Cell>();
      incidentCells[1] = new HashSet<Cell>();
      for (int i = 0; i < borderingCells.length; i++)
      {
         Cell bCell = borderingCells[i];
         if (!cells[cell.getDimension() - 1].contains(bCell))
            cells[cell.getDimension() - 1].add(bCell);

         incidentCells[0].add(bCell);

         if (temporaryIncidenceGraph.containsKey(bCell))
         {
            temporaryIncidenceGraph.get(bCell)[1].add(cell);
         }
         else
         {
            Set<Cell>[] bIncidentCells = (HashSet<Cell>[]) new HashSet[2];
            bIncidentCells[0] = new HashSet<Cell>();
            bIncidentCells[1] = new HashSet<Cell>();

            temporaryIncidenceGraph.put(bCell, bIncidentCells);
         }
      }
      if (!temporaryIncidenceGraph.containsValue(cell))
         temporaryIncidenceGraph.put(cell, incidentCells);
      else
      {
         temporaryIncidenceGraph.get(cell)[0].addAll(incidentCells[0]);
      }
      // incidence graph block===

      return cell;
   }

   @Override
   public Cell addCell(CellList borderingCells)
   {
      // TODO Auto-generated method stub
      Cell[] cellArr = new Cell[borderingCells.size()];
      borderingCells.toArray(cellArr);
      return addCell(cellArr);
   }

   @Override
   public boolean connectCells()
   {
      if (temporaryIncidenceGraph != null)
      {
         Set<Cell> newCells = temporaryIncidenceGraph.keySet();
         Iterator<Cell> it = newCells.iterator();

         ArrayList<Cell> nCells = new ArrayList<Cell>();
         while (it.hasNext())
         {
            Cell cell = it.next();
            if (cell.getDimension() == maxDimension)
               nCells.add(cell);
         }

         if (nCells.size() > 0)
         {
            Map<Cell, ArrayList<Dart>> cellDarts = new IdentityHashMap<Cell, ArrayList<Dart>>();
            //
            // if(this.darts.size()>0)
            // {
            // //TODO: prefill cellDarts with existing information to sew new stuff to old stuff
            // for(int i = 0; i <this.darts.size();i++)
            // {
            // Dart d = this.darts.get(i);
            // if(d.isOrphan()>-1)
            // {
            //
            // }
            // }
            // }

            ArrayList<Dart> newDartLst = new ArrayList<Dart>();

            for (int c = 0; c < nCells.size(); c += 1)
            {
               newDartLst.addAll(this.addNCell(nCells.get(c), temporaryIncidenceGraph, cellDarts));
            }

            boolean suc = check(newDartLst);
            if (suc)
               temporaryIncidenceGraph = null;
            return suc;
         }
         else
            return false;
      }
      else
         return false;

   }

   @Override
   public CellList getAllCells(int dimension)
   {
      return new CellList(cells[dimension]);
   }

   @Override
   public GMap intersect(HyperPlane plane)
   {
      GMap gmo = new GMap(maxDimension - 1);

      if (plane.normal.length == 3)
      {

         Map<Cell, Cell> oldNewCells = new IdentityHashMap<Cell, Cell>();
         Map<Cell, ArrayList<Dart>> orphanDartsAtCell = new IdentityHashMap<Cell, ArrayList<Dart>>();
         ArrayList<Dart> drtLst = new ArrayList<Dart>(darts);

         for (int k = 0; k < drtLst.size(); k++)
         {
            Dart d = drtLst.get(k);

            if (!d.isVisited())
            {
               Dart d2 = d.getInvolution(0);

               d.visit();
               if (d != d2)
               {
                  d2.visit();
                  Cell edgeCell = d.getIncidentCell(1);
                  Cell nVertexCell = null;
                  if (oldNewCells.containsKey(edgeCell))
                  {
                     nVertexCell = oldNewCells.get(edgeCell);
                  }
                  else
                  {
                     double[] p1 = d.getVertexCoordinates();
                     double[] p2 = d2.getVertexCoordinates();

                     double lambda = Base.cutPlaneEdgeV2(plane.normal, plane.position,
                           p1, p2);

                     if (lambda > -TOLERANCE_EPS && lambda < 1 + TOLERANCE_EPS)
                     {
                        double[] p = new double[p1.length];
                        for (int i = 0; i < p.length; i++)
                           p[i] = p1[i] + lambda * (p2[i] - p1[i]);

                        nVertexCell = new Cell(p);
                        oldNewCells.put(edgeCell, nVertexCell);
                        gmo.cells[0].add(nVertexCell);
                     }
                  }
                  if (nVertexCell != null)
                  {
                     // hier schnittpunkt
                     nVertexCell.setData(edgeCell.getData());
                     Dart newDart = new Dart(maxDimension - 1);
                     Cell[] incidentCells = newDart.getIncidentCells();
                     incidentCells[0] = nVertexCell;
                     // find incident cells
                     if (nVertexCell.getDart() == null)
                        nVertexCell.setDart(newDart);

                     if (orphanDartsAtCell.containsKey(nVertexCell))
                     {
                        orphanDartsAtCell.get(nVertexCell).add(newDart);
                     }
                     else
                     {
                        ArrayList<Dart> drts = new ArrayList<Dart>();
                        drts.add(newDart);
                        orphanDartsAtCell.put(nVertexCell, drts);
                     }

                     for (int i = 1; i < incidentCells.length; i++)
                     {
                        Cell iCell = d.getIncidentCell(i + 1);
                        Cell nICell = null;
                        if (oldNewCells.containsKey(iCell))
                        {
                           nICell = oldNewCells.get(iCell);
                        }
                        else
                        {
                           nICell = new Cell(i);
                           nICell.setDart(newDart);

                           gmo.cells[i].add(nICell);

                           oldNewCells.put(iCell, nICell);

                           nICell.setData(iCell.getData());
                        }

                        if (i == 1)
                        {
                           if (orphanDartsAtCell.containsKey(nICell))
                           {
                              orphanDartsAtCell.get(nICell).add(newDart);
                           }
                           else
                           {
                              ArrayList<Dart> drts = new ArrayList<Dart>();
                              drts.add(newDart);
                              orphanDartsAtCell.put(nICell, drts);
                           }
                        }

                        incidentCells[i] = nICell;
                     }

                     newDart.setIncidentCells(incidentCells);

                     newDart.setIncidentCells(incidentCells);

                     // here sewing: not time consuming
                     for (int i = 0; i <= newDart.getDimension(); i++)
                     {
                        ArrayList<Dart> drtToCheck = new ArrayList<Dart>();
                        if (i == 0)
                           drtToCheck = orphanDartsAtCell.get(newDart.getIncidentCell(1));
                        else
                           drtToCheck = orphanDartsAtCell.get(newDart.getIncidentCell(0));

                        // System.out.println(drtToCheck.size());

                        if (drtToCheck.size() > 0)
                        {
                           boolean found = false;
                           int nD = 0;

                           while (!found && nD < drtToCheck.size())
                           {
                              Dart drt2 = drtToCheck.get(nD);
                              if (newDart.isInvolution(drt2, i))
                              {
                                 newDart.sew(drt2, i);

                                 if (drt2.isOrphan() == -1)
                                 {
                                    drtToCheck.remove(drt2);
                                 }

                                 found = true;
                              }
                              else
                                 nD++;
                           }
                        }

                     }

                     int odn = newDart.isOrphan();
                     if (odn > -1)
                     {
                        boolean doAdd = false;
                        if (odn == 0)
                        {
                           if (orphanDartsAtCell.containsKey(newDart.getIncidentCell(1)))
                           {
                              orphanDartsAtCell.get(newDart.getIncidentCell(1)).add(newDart);
                           }
                           else
                           {
                              ArrayList<Dart> drts = new ArrayList<Dart>();
                              drts.add(newDart);
                              orphanDartsAtCell.put(newDart.getIncidentCell(1), drts);
                           }

                           for (int i = 1; i <= newDart.getDimension(); i++)
                           {
                              doAdd = newDart == newDart.getInvolution(i);
                           }
                        }

                        if (odn >= 1 || doAdd)
                        {
                           if (orphanDartsAtCell.containsKey(newDart.getIncidentCell(0)))
                           {
                              orphanDartsAtCell.get(newDart.getIncidentCell(0)).add(newDart);
                           }
                           else
                           {
                              ArrayList<Dart> drts = new ArrayList<Dart>();
                              drts.add(newDart);
                              orphanDartsAtCell.put(newDart.getIncidentCell(0), drts);
                           }
                        }
                     }

                     gmo.darts.add(newDart);

                  }

               }
            }
         }

         for (int k = 0; k < drtLst.size(); k++)
         {
            drtLst.get(k).unVisit();
         }
      }
      return gmo;
   }

   @Override
   public CellList CutByHyperPlane(HyperPlane plane) throws NotImplementedException
   {
      // TODO cut all cells and replace cutted cells by splitted cells
      throw new NotImplementedException("Not implemented yet");
   }

   @Override
   public GMap intersect(GMap gmap) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet");
   }

   @SuppressWarnings(
   { "unchecked" })
   @Override
   public GMap intersect(Cell ncell)
   {
      CellList cellVertices = ncell.getVertices();
      CellList mapVertices = this.getVertexCells();

      // connectivity function
      // association between on vertex v_i to the n-cell c_j incident to its nearest neighbor v_i is in c_j
      Map<Cell, Cell> connectivityFunction = new IdentityHashMap<Cell, Cell>();

      double[][] vertexArrCell = new double[cellVertices.size()][];
      for (int i = 0; i < cellVertices.size(); i++)
      {
         vertexArrCell[i] = cellVertices.get(i).getVertex();
         // System.out.println("cell vertex: " + Arrays.toString(vertexArrCell[i]));
      }

      double[][] vertexArrMap = new double[mapVertices.size()][];
      for (int i = 0; i < mapVertices.size(); i++)
      {
         vertexArrMap[i] = mapVertices.get(i).getVertex();
         // System.out.println("map vertex: " + Arrays.toString(vertexArrMap[i]));
      }

      // nearest neighbor
      NearestNeighbor nn_Map = new NearestNeighbor(vertexArrMap);
      NearestNeighbor nn_Cell = new NearestNeighbor(vertexArrCell);

      if (maxDimension == 2 && ncell.getDimension() == 2)
      {
         // 2D case intersection by a polygon
         // System.out.println(cellVertices.size());
         // TODO: connectivity function
         for (int i = 0; i < cellVertices.size(); i++)
         {
            double[] cellVrtxCoord = vertexArrCell[i];
            int id = nn_Map.getIndex(cellVrtxCoord);

            // System.out.println("NN-vertex: " + Arrays.toString(mapVertices.get(id).getVertex()) + " for "
            // + Arrays.toString(cellVrtxCoord));

            CellList incidentNCells = mapVertices.get(id).getNeighboringCells(2);
            int j = 0;
            while (j < incidentNCells.size() && !incidentNCells.get(j).containsPoint(cellVrtxCoord))
            {
               j++;
            }

            if (j < incidentNCells.size())
            {
               // System.out.println("incident case " + j);
               connectivityFunction.put(cellVertices.get(i), incidentNCells.get(j));
            }
            else
            {

               // TODO: point is not inside the cells incident to the NN
               // now test all cells along the ray between NN and test point
               Cell NNvertex = mapVertices.get(id);
               double[] nnCoords = NNvertex.getVertex();
               double[] rayDir = new double[cellVrtxCoord.length];
               for (int k = 0; k < cellVrtxCoord.length; k++)
                  rayDir[k] = cellVrtxCoord[k] - nnCoords[k];

               Cell nextCell = null;
               Cell oldCell = null;
               for (int k = 0; k < incidentNCells.size(); k++)
               {
                  Cell curFacet = incidentNCells.get(k);
                  oldCell = curFacet;
                  CellList fEdges = curFacet.getBorderingCells();
                  for (int l = 0; l < fEdges.size(); l++)
                  {
                     Cell edge = fEdges.get(l);
                     if (edge.getDart().getVertex() != NNvertex
                           && edge.getDart().getInvolution(0).getVertex() != NNvertex)
                     {

                        double[] egdeDir = new double[edge.getDart().getVertex().getVertex().length];
                        for (int m = 0; m < egdeDir.length; m++)
                           egdeDir[m] = edge.getDart().getInvolution(0).getVertex().getVertex()[m] -
                                 edge.getDart().getVertex().getVertex()[m];

                        // System.out.println(Arrays.toString(nnCoords)+"->"+Arrays.toString(rayDir)+" X
                        // "+Arrays.toString(edge.getDart().getVertex().getVertex())+"->"+Arrays.toString(egdeDir));
                        double[] lambda = Base.cutLineLine(nnCoords, rayDir,
                              edge.getDart().getVertex().getVertex(), egdeDir,
                              GMap.TOLERANCE_EPS);
                        // double[] lambda = Base.cutLineLine_C(nnCoords, rayDir,
                        // edge.getDart().getVertex().getVertex(), egdeDir);

                        // System.out.println(">>lambda: " + Arrays.toString(lambda));
                        if (!Double.isNaN(lambda[0]) && lambda[0] > 0 && lambda[0] < 1 &&
                              lambda[1] > GMap.TOLERANCE_EPS && lambda[1] < 1 + GMap.TOLERANCE_EPS)
                        {
                           if (edge.getDart() != edge.getDart().getInvolution(2))
                           {
                              if (edge.getDart().getIncidentCell(2) == curFacet)
                              {
                                 nextCell = edge.getDart().getInvolution(2).getIncidentCell(2);
                              }
                              else
                              {
                                 nextCell = edge.getDart().getIncidentCell(2);
                              }
                              break;
                           }
                        }
                     }
                  }
                  if (nextCell != null)
                     break;
               }

               while (nextCell != null && !nextCell.containsPoint(cellVrtxCoord))
               {
                  // System.out.println(">>traversing cell: "+nextCell + "(old cell: "+oldCell+")");
                  CellList fEdges = nextCell.getBorderingCells();
                  Cell curCell = nextCell;
                  nextCell = null;
                  for (int l = 0; l < fEdges.size(); l++)
                  {
                     Cell edge = fEdges.get(l);
                     if (edge.getDart().getIncidentCell(2) != oldCell
                           && edge.getDart().getInvolution(2).getIncidentCell(2) != oldCell)
                     {
                        // System.out.println(Arrays.toString()+"<>"+Arrays.toString()+" X
                        // "+Arrays.toString()+"<>"+Arrays.toString());
                        double[] egdeDir = new double[edge.getDart().getVertex().getVertex().length];
                        for (int m = 0; m < egdeDir.length; m++)
                           egdeDir[m] = edge.getDart().getInvolution(0).getVertex().getVertex()[m] -
                                 edge.getDart().getVertex().getVertex()[m];
                        double[] lambda = Base.cutLineLine(nnCoords, rayDir,
                              edge.getDart().getVertex().getVertex(), egdeDir,
                              GMap.TOLERANCE_EPS);
                        // double[] lambda = Base.cutLineLine_C(nnCoords, rayDir,
                        // edge.getDart().getVertex().getVertex(), egdeDir);

                        // System.out.println(">>lambda: " + Arrays.toString(lambda));
                        if (!Double.isNaN(lambda[0]) && lambda[0] > 0 && lambda[0] < 1 &&
                              lambda[1] > GMap.TOLERANCE_EPS && lambda[1] < 1 + GMap.TOLERANCE_EPS)
                        {
                           if (edge.getDart() != edge.getDart().getInvolution(2))
                           {
                              // System.out.println(">>cutting point found)");
                              if (edge.getDart().getIncidentCell(2) == curCell)
                              {
                                 nextCell = edge.getDart().getInvolution(2).getIncidentCell(2);
                              }
                              else
                              {
                                 nextCell = edge.getDart().getIncidentCell(2);
                              }
                              break;
                           }
                        }
                     }
                  }
                  oldCell = curCell;
               }

               if (nextCell == null)
               {
                  // System.out.println("No case");
                  connectivityFunction.put(cellVertices.get(i), NNvertex);
               }
               else
               {
                  // System.out.println("Other case");
                  connectivityFunction.put(cellVertices.get(i), nextCell);
               }
            }
         }

         for (int i = 0; i < mapVertices.size(); i++)
         {
            double[] cellVrtxCoord = vertexArrMap[i];

            Cell vc = mapVertices.get(i);

            if (ncell.containsPoint(cellVrtxCoord))
            {
               if (vc != null)
               {
                  connectivityFunction.put(vc, ncell);
               }

            }
            else
            {
               int id = nn_Cell.getIndex(cellVrtxCoord);
               connectivityFunction.put(vc, cellVertices.get(id));
            }
         }

         // print connectivity function
         // Set<Cell> keys = connectivityFunction.keySet();
         // Iterator<Cell> it = keys.iterator();
         // System.out.println("connectivity function:");
         // while (it.hasNext())
         // {
         // Cell key = it.next();
         // Cell value = connectivityFunction.get(key);
         // CellList valueVerts = value.getVertices();
         // System.out.println(Arrays.toString(key.getVertex()) + " in :");
         // for (int i = 0; i < valueVerts.size(); i++)
         // System.out.println("\t" + Arrays.toString(valueVerts.get(i).getVertex()));
         // }

         // TODO: 2D intersection - create post split list
         GMap gmo_return = new GMap(2);

         // test, whether all vertices of the cell are in the same cell

         int vId = 0;
         Cell inCell = connectivityFunction.get(cellVertices.get(vId));
         boolean allIn = inCell.getDimension() == 2;
         while (allIn && vId < cellVertices.size() - 1)
         {
            allIn = inCell == connectivityFunction.get(cellVertices.get(vId++));
         }

         if (allIn)
         {
            System.out.println("all in case");
            Map<Cell, Cell> oldNewCells = new IdentityHashMap<Cell, Cell>();
            Map<Cell, Set<Cell>[]> newIGraph = new IdentityHashMap<Cell, Set<Cell>[]>();

            Cell newFacet = new Cell(2);
            newFacet.setData(inCell.getData());

            oldNewCells.put(ncell, newFacet);

            Set<Cell>[] incidentCellsF = (HashSet<Cell>[]) new HashSet[2];
            incidentCellsF[0] = new HashSet<Cell>();
            incidentCellsF[1] = new HashSet<Cell>();

            newIGraph.put(newFacet, incidentCellsF);

            CellList bEdges = ncell.getBorderingCells();
            for (int e = 0; e < bEdges.size(); e++)
            {
               Cell oldEdge = bEdges.get(e);

               // System.out.println(e+" "+oldEdge+":
               // "+Arrays.toString(oldEdge.getDart().getVertexCoordinates())+"<>"+Arrays.toString(oldEdge.getDart().getInvolution(0).getVertexCoordinates()));

               Cell newEdge = oldNewCells.get(oldEdge);
               Set<Cell>[] incidentCellsE;

               if (newEdge == null)
               {
                  newEdge = new Cell(1);
                  oldNewCells.put(oldEdge, newEdge);

                  incidentCellsE = (HashSet<Cell>[]) new HashSet[2];
                  incidentCellsE[0] = new HashSet<Cell>();
                  incidentCellsE[1] = new HashSet<Cell>();

                  incidentCellsE[1].add(newFacet);

                  newIGraph.put(newEdge, incidentCellsE);

                  Cell oldEVC1 = oldEdge.getDart().getVertex();
                  Cell oldEVC2 = oldEdge.getDart().getInvolution(0).getVertex();

                  Cell newEVC1 = oldNewCells.get(oldEVC1);
                  Set<Cell>[] incidentCellsV1;

                  if (newEVC1 == null)
                  {
                     newEVC1 = new Cell(oldEVC1.getVertex());
                     oldNewCells.put(oldEVC1, newEVC1);

                     incidentCellsV1 = (HashSet<Cell>[]) new HashSet[2];
                     incidentCellsV1[0] = new HashSet<Cell>();
                     incidentCellsV1[1] = new HashSet<Cell>();

                     incidentCellsV1[1].add(newEdge);

                     newIGraph.put(newEVC1, incidentCellsV1);
                  }
                  else
                  {

                     incidentCellsV1 = newIGraph.get(newEVC1);
                     if (incidentCellsV1 == null)
                     {
                        incidentCellsV1 = (HashSet<Cell>[]) new HashSet[2];
                        incidentCellsV1[0] = new HashSet<Cell>();
                        incidentCellsV1[1] = new HashSet<Cell>();
                        newIGraph.put(newEVC1, incidentCellsV1);
                     }
                     incidentCellsV1[1].add(newEdge);
                  }

                  Cell newEVC2 = oldNewCells.get(oldEVC2);
                  Set<Cell>[] incidentCellsV2;

                  if (newEVC2 == null)
                  {
                     newEVC2 = new Cell(oldEVC2.getVertex());
                     oldNewCells.put(oldEVC2, newEVC2);

                     incidentCellsV2 = (HashSet<Cell>[]) new HashSet[2];
                     incidentCellsV2[0] = new HashSet<Cell>();
                     incidentCellsV2[1] = new HashSet<Cell>();

                     incidentCellsV2[1].add(newEdge);

                     newIGraph.put(newEVC2, incidentCellsV2);
                  }
                  else
                  {
                     incidentCellsV2 = newIGraph.get(newEVC2);
                     if (incidentCellsV2 == null)
                     {
                        incidentCellsV2 = (HashSet<Cell>[]) new HashSet[2];
                        incidentCellsV2[0] = new HashSet<Cell>();
                        incidentCellsV2[1] = new HashSet<Cell>();
                        newIGraph.put(newEVC2, incidentCellsV1);
                     }
                     incidentCellsV2[1].add(newEdge);
                  }

                  incidentCellsE[0].add(newEVC1);
                  incidentCellsE[0].add(newEVC2);
               }
               else
               {
                  incidentCellsE = newIGraph.get(newEdge);
                  incidentCellsE[1].add(newFacet);
               }

               incidentCellsF[0].add(newEdge);
            }
            gmo_return.buildGMapByIncidenceGraphGeneralized(newIGraph);
         }
         else
         {
            Map<Cell, ArrayList<Cell[]>> postSplitList = new IdentityHashMap<Cell, ArrayList<Cell[]>>();
            Map<Cell, Cell> oldNewCells = new IdentityHashMap<Cell, Cell>();
            CellList edges = ncell.getBorderingCells();
            for (int e = 0; e < edges.size(); e++)
            {
               Cell u = edges.get(e).getDart().getVertex();
               Cell v = edges.get(e).getDart().getInvolution(0).getVertex();

               // System.out.println("Edge " + Arrays.toString(u.getVertex()) + "->" + Arrays.toString(v.getVertex()));

               // System.out.println("Cell around u: " + connectivityFunction.get(u));
               // System.out.println("Cell around v: " + connectivityFunction.get(v));

               Cell[] newEdgeVertices = new Cell[2];
               if (connectivityFunction.get(u).getDimension() == maxDimension
                     && connectivityFunction.get(v).getDimension() == maxDimension
                     && connectivityFunction.get(u) == connectivityFunction.get(v))
               {
                  System.out.println("this only happens when both in the same cell");
                  // Warning: This approach assumes that when both points of an edge are inside the same polygon
                  // , this edge never intersects the edges of the polygon. This is only true for "good shaped"
                  // (at least convex) polygons

                  if (oldNewCells.containsKey(u))
                  {
                     newEdgeVertices[0] = oldNewCells.get(u);
                  }
                  else
                  {
                     newEdgeVertices[0] = new Cell(u.getVertex());
                     oldNewCells.put(u, newEdgeVertices[0]);
                  }
                  if (oldNewCells.containsKey(v))
                  {
                     newEdgeVertices[1] = oldNewCells.get(v);
                  }
                  else
                  {
                     newEdgeVertices[1] = new Cell(v.getVertex());
                     oldNewCells.put(v, newEdgeVertices[1]);
                  }

                  Cell cell = connectivityFunction.get(u);

                  // System.out.println(cell+": "+Arrays.toString(newEdgeVertices));
                  if (postSplitList.containsKey(cell))
                  {
                     postSplitList.get(cell).add(newEdgeVertices);
                  }
                  else
                  {
                     ArrayList<Cell[]> lst = new ArrayList<Cell[]>();
                     lst.add(newEdgeVertices);
                     postSplitList.put(cell, lst);
                  }

               }
               else
                  if (connectivityFunction.get(u).getDimension() == 0
                        || connectivityFunction.get(v).getDimension() == 0)
                  {
                     System.out.println("out case");

                     Cell curCell = null;
                     Cell edge = null;
                     Cell startVrtx = null;
                     Cell endVrtx = null;

                     if (connectivityFunction.get(u).getDimension() != 0)
                     {
                        curCell = connectivityFunction.get(u);
                        if (oldNewCells.containsKey(u))
                        {
                           startVrtx = oldNewCells.get(u);
                        }
                        else
                        {
                           startVrtx = new Cell(u.getVertex());
                           oldNewCells.put(u, startVrtx);
                        }
                     }
                     else
                        if (connectivityFunction.get(v).getDimension() != 0)
                        {
                           curCell = connectivityFunction.get(v);
                           if (oldNewCells.containsKey(v))
                           {
                              startVrtx = oldNewCells.get(v);
                           }
                           else
                           {
                              startVrtx = new Cell(v.getVertex());
                              oldNewCells.put(v, startVrtx);
                           }
                        }

                     if (startVrtx == null)
                     {

                        Dart oldDart = connectivityFunction.get(u).getDart();
                        Dart d1 = connectivityFunction.get(u).getDart();
                        double[] w = new double[d1.getVertex().getVertex().length];
                        while (d1 != d1.getInvolution(2))
                        {
                           d1 = d1.getInvolution(2).getInvolution(1);
                        }
                        edge = d1.getIncidentCell(1);
                        curCell = d1.getIncidentCell(2);
                        // System.out.println(">>Edge: "+
                        // Arrays.toString(d1.getVertexCoordinates())+"->"+Arrays.toString(d1.getInvolution(0).getVertexCoordinates()));
                        w = Base.cutEdgeEdge(u.getVertex(), v.getVertex(),
                              d1.getVertexCoordinates(),
                              d1.getInvolution(0).getVertexCoordinates(),
                              GMap.TOLERANCE_EPS);

                        if (w == null)
                        {
                           Dart d2 = d1.getInvolution(1);
                           while (d2 != d2.getInvolution(2))
                           {
                              d2 = d2.getInvolution(2).getInvolution(1);
                           }
                           edge = d2.getIncidentCell(1);
                           curCell = d2.getIncidentCell(2);
                           // System.out.println(">>Edge: "+
                           // Arrays.toString(d2.getVertexCoordinates())+"->"+Arrays.toString(d2.getInvolution(0).getVertexCoordinates()));
                           w = Base.cutEdgeEdge(u.getVertex(), v.getVertex(),
                                 d2.getVertexCoordinates(),
                                 d2.getInvolution(0).getVertexCoordinates(),
                                 GMap.TOLERANCE_EPS);

                           while (w == null && d1.getInvolution(0) != oldDart && d2.getInvolution(0) != oldDart)
                           {
                              if (d1.getInvolution(0) != oldDart)
                              {
                                 d1 = d1.getInvolution(0).getInvolution(1);
                                 while (d1 != d1.getInvolution(2))
                                 {
                                    d1 = d1.getInvolution(2).getInvolution(1);
                                 }

                                 edge = d1.getIncidentCell(1);
                                 curCell = d1.getIncidentCell(2);
                                 // System.out.println(">>Edge: "+
                                 // Arrays.toString(d1.getVertexCoordinates())+"->"+Arrays.toString(d1.getInvolution(0).getVertexCoordinates()));
                                 w = Base.cutEdgeEdge(u.getVertex(), v.getVertex(),
                                       d1.getVertexCoordinates(),
                                       d1.getInvolution(0).getVertexCoordinates(),
                                       GMap.TOLERANCE_EPS);
                              }
                              if (w == null)
                              {
                                 if (d2.getInvolution(0) != oldDart)
                                 {
                                    d2 = d2.getInvolution(0).getInvolution(1);
                                    while (d2 != d2.getInvolution(2))
                                    {
                                       d2 = d2.getInvolution(2).getInvolution(1);
                                    }
                                    edge = d2.getIncidentCell(1);
                                    curCell = d2.getIncidentCell(2);
                                    // System.out.println(">>Edge: "+
                                    // Arrays.toString(d2.getVertexCoordinates())+"->"+Arrays.toString(d2.getInvolution(0).getVertexCoordinates()));
                                    w = Base.cutEdgeEdge(u.getVertex(), v.getVertex(),
                                          d2.getVertexCoordinates(),
                                          d2.getInvolution(0).getVertexCoordinates(),
                                          GMap.TOLERANCE_EPS);
                                 }
                              }
                           }
                        }
                        if (w != null)
                           startVrtx = new Cell(w);
                        // System.out.println(">>>>"+startVrtx.print());
                     }

                     if (startVrtx != null && curCell != null)
                     {
                        CellList cellEdges = curCell.getBorderingCells();
                        Cell cutEdge = edge;
                        while (curCell != connectivityFunction.get(v));
                        {
                           if (cutEdge != null)
                              cellEdges.remove(cutEdge);
                        }

                        double[] w = null;
                        for (int e1 = 0; e1 < cellEdges.size(); e1++)
                        {
                           Cell cedge = cellEdges.get(e1);

                           w = Base.cutEdgeEdge(u.getVertex(), v.getVertex(),
                                 cedge.getDart().getVertexCoordinates(),
                                 cedge.getDart().getInvolution(0).getVertexCoordinates(),
                                 GMap.TOLERANCE_EPS);

                           if (w != null)
                           {
                              endVrtx = new Cell(w);
                              cutEdge = cedge;
                              break;
                           }
                        }

                        if (cutEdge != null && endVrtx != null)
                        {
                           System.out.println(">> with intersections - not supported yet");
                           // newEdgeVertices[0] = startVrtx;
                           // newEdgeVertices[1] = endVrtx;
                           //
                           // if (postSplitList.containsKey(curCell))
                           // {
                           // postSplitList.get(curCell).add(newEdgeVertices);
                           // }
                           // else
                           // {
                           // ArrayList<Cell[]> lst = new ArrayList<Cell[]>();
                           // lst.add(newEdgeVertices);
                           // postSplitList.put(curCell, lst);
                           // }
                           //
                           // System.out.println("cut edge " + Arrays.toString(newEdgeVertices[0].getVertex()) + "->"
                           // + Arrays.toString(newEdgeVertices[1].getVertex()));
                           //
                           // Cell[] innerEdge = null;
                           //
                           // // add new segment for inner edge
                           // if (!connectivityFunction.containsKey(cutEdge.getDart().getVertex()) &&
                           // connectivityFunction.containsKey(cutEdge.getDart().getInvolution(0).getVertex()))
                           // {
                           // System.out.println("edge ends in cell");
                           // Cell tmp = cutEdge.getDart().getInvolution(0).getVertex();
                           // Cell iV;
                           // if (oldNewCells.containsKey(tmp))
                           // {
                           // iV = oldNewCells.get(tmp);
                           // }
                           // else
                           // {
                           // iV = new Cell(tmp.getVertex());
                           // oldNewCells.put(tmp, iV);
                           // }
                           // innerEdge = new Cell[2];
                           // innerEdge[0] = endVrtx;
                           // innerEdge[1] = iV;
                           //
                           // postSplitList.get(curCell).add(innerEdge);
                           //
                           // }
                           // else
                           // if (connectivityFunction.containsKey(cutEdge.getDart().getVertex()) &&
                           // !connectivityFunction.containsKey(cutEdge.getDart().getInvolution(0).getVertex()))
                           // {
                           // System.out.println("edge starts in cell");
                           // Cell tmp = cutEdge.getDart().getVertex();
                           // Cell iV;
                           // if (oldNewCells.containsKey(tmp))
                           // {
                           // iV = oldNewCells.get(tmp);
                           // }
                           // else
                           // {
                           // iV = new Cell(tmp.getVertex());
                           // oldNewCells.put(tmp, iV);
                           // }
                           // innerEdge = new Cell[2];
                           // innerEdge[0] = endVrtx;
                           // innerEdge[1] = iV;
                           // postSplitList.get(curCell).add(innerEdge);
                           // }
                           // else
                           // if (!connectivityFunction.containsKey(cutEdge.getDart().getVertex()) &&
                           // !connectivityFunction.containsKey(cutEdge.getDart().getInvolution(0).getVertex()))
                           // {
                           // // none is in, but cut ...
                           // System.out.println("edge cuts cell");
                           // }
                           //
                           // startVrtx = endVrtx;
                           //
                           // if (cutEdge.getDart().getIncidentCell(2) != curCell)
                           // {
                           // curCell = cutEdge.getDart().getIncidentCell(2);
                           // }
                           // else
                           // {
                           // curCell = cutEdge.getDart().getInvolution(2).getIncidentCell(2);
                           // }
                           //
                           // startVrtx = endVrtx;
                           //
                           // if (cutEdge.getDart().getIncidentCell(2) != curCell)
                           // {
                           // curCell = cutEdge.getDart().getIncidentCell(2);
                           // }
                           // else
                           // {
                           // curCell = cutEdge.getDart().getInvolution(2).getIncidentCell(2);
                           // }
                           //
                           // if (innerEdge != null)
                           // {
                           // System.out.println("inner edge " + Arrays.toString(innerEdge[0].getVertex()) + "->"
                           // + Arrays.toString(innerEdge[1].getVertex()));
                           //
                           // if (postSplitList.containsKey(curCell))
                           // {
                           // ArrayList<Cell[]> lst = postSplitList.get(curCell);
                           // lst.add(newEdgeVertices);
                           // lst.add(innerEdge);
                           // }
                           // else
                           // {
                           // ArrayList<Cell[]> lst = new ArrayList<Cell[]>();
                           // lst.add(newEdgeVertices);
                           // lst.add(innerEdge);
                           // postSplitList.put(curCell, lst);
                           // }
                           // }
                           // // else
                           // // System.out.println("keine inner edge?");
                           //
                        }

                     }
                     else
                     {
                        // no cuts, this edge lies completely outside, nothing to do here
                     }
                  }
                  else
                  {
                     System.out.println("general case");

                     Cell curCell = connectivityFunction.get(u);

                     Cell startVrtx = null;
                     if (oldNewCells.containsKey(u))
                     {
                        startVrtx = oldNewCells.get(u);
                     }
                     else
                     {
                        startVrtx = new Cell(u.getVertex());
                        oldNewCells.put(u, startVrtx);
                     }
                     // System.out.println(">> find start " + curCell);
                     // System.out.println(">> " + connectivityFunction.get(v));
                     // System.out.println(">> " + (curCell != connectivityFunction.get(v)));

                     Cell endVrtx = null;
                     Cell cutEdge = null;

                     while (curCell != connectivityFunction.get(v))
                     {
                        // System.out.println("curCell != connectivityFunction.get(v)");
                        CellList cellEdges = curCell.getBorderingCells();
                        if (cutEdge != null)
                           cellEdges.remove(cutEdge);

                        double[] w = null;
                        for (int e1 = 0; e1 < cellEdges.size(); e1++)
                        {
                           Cell cedge = cellEdges.get(e1);

                           w = Base.cutEdgeEdge(u.getVertex(), v.getVertex(),
                                 cedge.getDart().getVertexCoordinates(),
                                 cedge.getDart().getInvolution(0).getVertexCoordinates(),
                                 GMap.TOLERANCE_EPS);

                           if (w != null)
                           {
                              endVrtx = new Cell(w);
                              cutEdge = cedge;
                              break;
                           }
                        }

                        if (cutEdge != null && endVrtx != null)
                        {
                           newEdgeVertices = new Cell[2];
                           newEdgeVertices[0] = startVrtx;
                           newEdgeVertices[1] = endVrtx;

                           // System.out
                           // .println("cut edge " + Arrays.toString(cutEdge.getDart().getVertexCoordinates()) + "->"
                           // + Arrays.toString(cutEdge.getDart().getInvolution(0).getVertexCoordinates()));
                           //
                           // System.out.println("new edge " + Arrays.toString(newEdgeVertices[0].getVertex()) + "->"
                           // + Arrays.toString(newEdgeVertices[1].getVertex()));

                           if (postSplitList.containsKey(curCell))
                           {
                              // System.out.println("-- cut edge" + curCell);
                              postSplitList.get(curCell).add(newEdgeVertices);

                              // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                              // for (int i = 0; i < segments.size(); i++)
                              // {
                              // System.out.println("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                              // + Arrays.toString(segments.get(i)[1].getVertex()));
                              // }
                           }
                           else
                           {
                              // System.out.println("++ cut edge" + curCell);
                              ArrayList<Cell[]> lst = new ArrayList<Cell[]>();
                              lst.add(newEdgeVertices);
                              postSplitList.put(curCell, lst);

                              // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                              // for (int i = 0; i < segments.size(); i++)
                              // {
                              // System.out.println("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                              // + Arrays.toString(segments.get(i)[1].getVertex()));
                              // }
                           }

                           Cell[] innerEdge = null;

                           // add new segment for inner edge
                           // System.out.println("Edge " + Arrays.toString(cutEdge.getDart().getVertex().getVertex()) +
                           // "->" + Arrays.toString(cutEdge.getDart().getInvolution(0).getVertex().getVertex()));
                           // System.out.println(connectivityFunction.get(cutEdge.getDart().getVertex()).getDimension()
                           // +"/"+connectivityFunction.get(cutEdge.getDart().getInvolution(0).getVertex()).getDimension());
                           if (connectivityFunction.get(cutEdge.getDart().getVertex()).getDimension() == 0 &&
                                 connectivityFunction.get(cutEdge.getDart().getInvolution(0).getVertex())
                                       .getDimension() != 0)
                           {
                              Cell tmp = cutEdge.getDart().getInvolution(0).getVertex();
                              // Cell neighbor = cutEdge.getDart().getInvolution(2).getIncidentCell(2);

                              Cell iV;
                              if (oldNewCells.containsKey(tmp))
                              {
                                 iV = oldNewCells.get(tmp);
                              }
                              else
                              {
                                 iV = new Cell(tmp.getVertex());
                                 oldNewCells.put(tmp, iV);
                              }

                              innerEdge = new Cell[3];
                              innerEdge[0] = endVrtx;
                              innerEdge[1] = iV;
                              innerEdge[2] = cutEdge;

                              // System.out.println("-- inner edge 1 " + curCell);
                              postSplitList.get(curCell).add(innerEdge);
                              // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                              // for (int i = 0; i < segments.size(); i++)
                              // {
                              // System.out.println("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                              // + Arrays.toString(segments.get(i)[1].getVertex()));
                              // }

                           }
                           else
                              if (connectivityFunction.get(cutEdge.getDart().getVertex()).getDimension() != 0 &&
                                    connectivityFunction.get(cutEdge.getDart().getInvolution(0).getVertex())
                                          .getDimension() == 0)
                              {
                                 Cell tmp = cutEdge.getDart().getVertex();
                                 // Cell neighbor = cutEdge.getDart().getInvolution(2).getIncidentCell(2);

                                 Cell iV;
                                 if (oldNewCells.containsKey(tmp))
                                 {
                                    iV = oldNewCells.get(tmp);
                                 }
                                 else
                                 {
                                    iV = new Cell(tmp.getVertex());
                                    oldNewCells.put(tmp, iV);
                                 }

                                 innerEdge = new Cell[3];
                                 innerEdge[0] = iV;
                                 innerEdge[1] = endVrtx;
                                 innerEdge[2] = cutEdge;

                                 // System.out.println("-- inner edge 2 " + curCell);
                                 postSplitList.get(curCell).add(innerEdge);
                                 // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                                 // for (int i = 0; i < segments.size(); i++)
                                 // {
                                 // System.out.println("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                                 // + Arrays.toString(segments.get(i)[1].getVertex()));
                                 // }

                              }
                              else
                              {
                                 // todo:
                                 // System.out.println("Hier wirds problematisch");
                                 innerEdge = new Cell[3];
                                 innerEdge[0] = endVrtx;
                                 innerEdge[1] = endVrtx;
                                 innerEdge[2] = cutEdge;
                                 postSplitList.get(curCell).add(innerEdge);
                              }
                           startVrtx = endVrtx;
                           // System.out.println("old curCell: "+ curCell);
                           if (cutEdge.getDart().getIncidentCell(2) != curCell)
                           {
                              curCell = cutEdge.getDart().getIncidentCell(2);
                           }
                           else
                           {
                              curCell = cutEdge.getDart().getInvolution(2).getIncidentCell(2);
                           }

                           // System.out.println("new curCell: "+ curCell);
                           if (innerEdge != null)
                           {
                              // System.out.println("inner edge " + Arrays.toString(innerEdge[0].getVertex()) + "->"
                              // + Arrays.toString(innerEdge[1].getVertex()));

                              if (postSplitList.containsKey(curCell))
                              {

                                 ArrayList<Cell[]> lst = postSplitList.get(curCell);
                                 lst.add(innerEdge);
                                 // System.out.println(">>neighbor exists " + lst.size());
                              }
                              else
                              {

                                 ArrayList<Cell[]> lst = new ArrayList<Cell[]>();
                                 lst.add(innerEdge);
                                 postSplitList.put(curCell, lst);
                                 // System.out.println(">>neighbor new " + lst.size());
                              }
                              // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                              // System.out.println("-- " + curCell + " " + segments.size());
                              // for (int i = 0; i < segments.size(); i++)
                              // {
                              // System.out
                              // .println("\t" + i + ":" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                              // + Arrays.toString(segments.get(i)[1].getVertex()));
                              // }

                           }

                        }
                        // System.out.println(">>> curCell:"+ curCell);
                     }

                     // System.out.println(">> while finished");

                     // last line segment in connectivityFunction.get(v)
                     if (oldNewCells.containsKey(v))
                     {
                        endVrtx = oldNewCells.get(v);
                     }
                     else
                     {
                        endVrtx = new Cell(v.getVertex());
                        oldNewCells.put(v, endVrtx);
                     }

                     newEdgeVertices = new Cell[2];
                     newEdgeVertices[0] = startVrtx;
                     newEdgeVertices[1] = endVrtx;

                     // System.out.println("new edge " + Arrays.toString(newEdgeVertices[0].getVertex()) + "->"
                     // + Arrays.toString(newEdgeVertices[1].getVertex()));
                     //
                     // System.out.println("new curCell: " + curCell);

                     if (postSplitList.containsKey(curCell))
                     {
                        // System.out.println("-- cut edge" + curCell);
                        postSplitList.get(curCell).add(newEdgeVertices);
                        // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                        // for (int i = 0; i < segments.size(); i++)
                        // {
                        // System.out.println("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                        // + Arrays.toString(segments.get(i)[1].getVertex()));
                        // }
                     }
                     else
                     {
                        // System.out.println("++ cut edge" + curCell);
                        ArrayList<Cell[]> lst = new ArrayList<Cell[]>();
                        lst.add(newEdgeVertices);
                        postSplitList.put(curCell, lst);

                        // ArrayList<Cell[]> segments = postSplitList.get(curCell);
                        // for (int i = 0; i < segments.size(); i++)
                        // {
                        // System.out.println("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                        // + Arrays.toString(segments.get(i)[1].getVertex()));
                        // }
                     }

                  }
            }

            // splitting
            System.out.println("<<<Splitting>>>");
            Set<Cell> splittedCells = postSplitList.keySet();

            Set<Cell> coveredCells = new HashSet<Cell>();;
            if (splittedCells.size() > 0)
            {
               Iterator<Cell> itSC = splittedCells.iterator();
               while (itSC.hasNext())
               {
                  Cell splittedCell = itSC.next();

                  ArrayList<Cell[]> segments = postSplitList.get(splittedCell);

                  ArrayList<Cell> startEndeEdgeL = new ArrayList<Cell>();
                  // Cell[] startEndeEdge = new Cell[2];
                  // System.out.println(splittedCell + ":");
                  // int j = 0;
                  for (int i = 0; i < segments.size(); i++)
                  {
                     // System.out.print("\t" + Arrays.toString(segments.get(i)[0].getVertex()) + "<->"
                     // + Arrays.toString(segments.get(i)[1].getVertex()));

                     if (segments.get(i).length == 3)
                     {
                        startEndeEdgeL.add(segments.get(i)[2]);
                        // startEndeEdge[j] = segments.get(i)[2];
                        // j++;
                        // System.out.println("\t\t replacing edge " + segments.get(i)[2]);
                     }
                     // else
                     // System.out.println("\t\t new edge");
                  }

                  if (startEndeEdgeL.size() == 2)
                  {
                     // System.out.println(">>general case: " + startEndeEdgeL.get(0) + " <> " + startEndeEdgeL.get(1));
                     // old and standart case

                     if (startEndeEdgeL.get(0) == startEndeEdgeL.get(1))
                     {
                        Cell[] newEdge = new Cell[2];
                        Cell oldEdge = startEndeEdgeL.get(0);
                        int k = 0;
                        for (int i = 0; i < segments.size(); i++)
                        {
                           if (segments.get(i).length == 3 && segments.get(i)[2] == startEndeEdgeL.get(0))
                           {
                              newEdge[k] = segments.get(i)[0];
                              k++;
                           }
                           // else
                           // if (segments.get(i).length == 3 && segments.get(i)[2] == startEndeEdgeL.get(1))
                           // {
                           // newEdge[1] = segments.get(i)[0];
                           // }
                        }

                        // System.out.println(">>>>"+newEdge[0]+"<>"+newEdge[1]);

                        Dart drt = oldEdge.getDart();
                        if (!splittedCells.contains(drt.getIncidentCell(2))
                              && !coveredCells.contains(drt.getIncidentCell(2)))
                        {
                           coveredCells.add(drt.getIncidentCell(2));
                        }

                        if (!splittedCells.contains(drt.getInvolution(2).getIncidentCell(2))
                              && !coveredCells.contains(drt.getInvolution(2).getIncidentCell(2)))
                        {
                           coveredCells.add(drt.getInvolution(2).getIncidentCell(2));
                        }
                        postSplitList.get(splittedCell).add(newEdge);
                     }
                     else
                     {
                        Dart eDart = startEndeEdgeL.get(0).getDart();
                        if (connectivityFunction.get(eDart.getVertex()).getDimension() == 0)
                        {
                           eDart = eDart.getInvolution(0);
                        }
                        if (eDart.getIncidentCell(2) != splittedCell)
                        {
                           eDart = eDart.getInvolution(2);
                        }

                        eDart = eDart.getInvolution(1);

                        Cell nextEdge = eDart.getIncidentCell(1);

                        while (nextEdge != startEndeEdgeL.get(1))
                        {
                           Cell[] newEdge = new Cell[2];

                           if (!splittedCells.contains(eDart.getInvolution(2).getIncidentCell(2)))
                              coveredCells.add(eDart.getInvolution(2).getIncidentCell(2));

                           if (oldNewCells.containsKey(eDart.getVertex()))
                           {
                              newEdge[0] = oldNewCells.get(eDart.getVertex());
                           }
                           else
                           {
                              newEdge[0] = new Cell(eDart.getVertex().getVertex());
                              oldNewCells.put(eDart.getVertex(), newEdge[0]);
                           }

                           if (oldNewCells.containsKey(eDart.getInvolution(0).getVertex()))
                           {
                              newEdge[1] = oldNewCells.get(eDart.getInvolution(0).getVertex());
                           }
                           else
                           {
                              newEdge[1] = new Cell(eDart.getInvolution(0).getVertex().getVertex());
                              oldNewCells.put(eDart.getInvolution(0).getVertex(), newEdge[1]);
                           }

                           postSplitList.get(splittedCell).add(newEdge);

                           eDart = eDart.getInvolution(0).getInvolution(1);
                           nextEdge = eDart.getIncidentCell(1);
                        }
                     }
                  }
                  else
                  {
                     // System.out.println(">>special case");
                     // special case
                     while (!startEndeEdgeL.isEmpty())
                     {
                        Cell startEdge = startEndeEdgeL.get(0);
                        Cell endEdge = null;

                        for (int i = 1; i < startEndeEdgeL.size(); i++)
                        {
                           if (startEndeEdgeL.get(i) == startEdge)
                           {
                              endEdge = startEndeEdgeL.get(i);
                              break;
                           }
                        }

                        if (endEdge != null)
                        {
                           // single edge case
                           Cell[] newEdge = new Cell[2];
                           Cell oldEdge = startEdge;
                           int k = 0;
                           for (int i = 0; i < segments.size(); i++)
                           {
                              if (segments.get(i).length == 3 && segments.get(i)[2] == startEndeEdgeL.get(0))
                              {
                                 newEdge[k] = segments.get(i)[0];
                                 k++;
                              }
                              // else
                              // if (segments.get(i).length == 3 && segments.get(i)[2] == startEndeEdgeL.get(1))
                              // {
                              // newEdge[1] = segments.get(i)[0];
                              // }
                           }

                           // System.out.println(">>>>"+newEdge[0]+"<>"+newEdge[1]);

                           Dart drt = oldEdge.getDart();
                           if (!splittedCells.contains(drt.getIncidentCell(2))
                                 && !coveredCells.contains(drt.getIncidentCell(2)))
                           {
                              coveredCells.add(drt.getIncidentCell(2));
                           }

                           if (!splittedCells.contains(drt.getInvolution(2).getIncidentCell(2))
                                 && !coveredCells.contains(drt.getInvolution(2).getIncidentCell(2)))
                           {
                              coveredCells.add(drt.getInvolution(2).getIncidentCell(2));
                           }

                           postSplitList.get(splittedCell).add(newEdge);
                        }
                        else
                        {
                           // inner edges case
                           Dart eDart = startEndeEdgeL.get(0).getDart();
                           if (connectivityFunction.get(eDart.getVertex()).getDimension() == 0)
                           {
                              eDart = eDart.getInvolution(0);
                           }
                           if (eDart.getIncidentCell(2) != splittedCell)
                           {
                              eDart = eDart.getInvolution(2);
                           }

                           eDart = eDart.getInvolution(1);

                           Cell nextEdge = eDart.getIncidentCell(1);
                           while (endEdge == null)
                           {

                              Cell[] newEdge = new Cell[2];
                              if (!splittedCells.contains(eDart.getInvolution(2).getIncidentCell(2)))
                                 coveredCells.add(eDart.getInvolution(2).getIncidentCell(2));

                              if (oldNewCells.containsKey(eDart.getVertex()))
                              {
                                 newEdge[0] = oldNewCells.get(eDart.getVertex());
                              }
                              else
                              {
                                 newEdge[0] = new Cell(eDart.getVertex().getVertex());
                                 oldNewCells.put(eDart.getVertex(), newEdge[0]);
                              }

                              if (oldNewCells.containsKey(eDart.getInvolution(0).getVertex()))
                              {
                                 newEdge[1] = oldNewCells.get(eDart.getInvolution(0).getVertex());
                              }
                              else
                              {
                                 newEdge[1] = new Cell(eDart.getInvolution(0).getVertex().getVertex());
                                 oldNewCells.put(eDart.getInvolution(0).getVertex(), newEdge[1]);
                              }

                              postSplitList.get(splittedCell).add(newEdge);

                              eDart = eDart.getInvolution(0).getInvolution(1);
                              nextEdge = eDart.getIncidentCell(1);

                              for (int i = 1; i < startEndeEdgeL.size(); i++)
                              {
                                 if (startEndeEdgeL.get(i) == nextEdge)
                                 {
                                    endEdge = startEndeEdgeL.get(i);
                                    break;
                                 }
                              }
                           }
                        }
                        startEndeEdgeL.remove(startEdge);
                        startEndeEdgeL.remove(endEdge);
                     }
                  }
               }

               // additional covered cells
               System.out.println("<<<Covering - partial>>>");
               Set<Cell> curCells = new HashSet<Cell>(coveredCells);
               while (!curCells.isEmpty())
               {
                  Iterator<Cell> itCC = curCells.iterator();
                  Set<Cell> tmp = new HashSet<Cell>();
                  while (itCC.hasNext())
                  {
                     CellList neighbors = itCC.next().getNeighboringCells();
                     for (int n = 0; n < neighbors.size(); n++)
                     {
                        Cell neighbor = neighbors.get(n);
                        if (!splittedCells.contains(neighbor) && !coveredCells.contains(neighbor))
                        {
                           coveredCells.add(neighbor);
                           tmp.add(neighbor);
                        }
                     }
                  }
                  curCells = new HashSet<Cell>(tmp);
               }

            }
            else
            {
               System.out.println("<<<Covering - global>>>");
               coveredCells = new HashSet<Cell>(this.getFacetCells().getList());
            }

            System.out.println("<<<building new GMap>>>");
            Map<Cell, Set<Cell>[]> newIGraph = new IdentityHashMap<Cell, Set<Cell>[]>();

            Iterator<Cell> itCC = coveredCells.iterator();
            System.out.println("<<<<covered cells");
            while (itCC.hasNext())
            {

               Cell oldFacet = itCC.next();

               Cell newFacet = new Cell(2);
               newFacet.setData(oldFacet.getData());
               oldNewCells.put(oldFacet, newFacet);

               Set<Cell>[] incidentCellsF = (HashSet<Cell>[]) new HashSet[2];
               incidentCellsF[0] = new HashSet<Cell>();
               incidentCellsF[1] = new HashSet<Cell>();

               newIGraph.put(newFacet, incidentCellsF);

               CellList bEdges = oldFacet.getBorderingCells();
               for (int e = 0; e < bEdges.size(); e++)
               {
                  Cell oldEdge = bEdges.get(e);

                  if (oldEdge.getDart().getVertex() != oldEdge.getDart().getInvolution(0).getVertex())
                  {
                     Cell newEdge = oldNewCells.get(oldEdge);
                     Set<Cell>[] incidentCellsE;

                     if (newEdge == null)
                     {
                        newEdge = new Cell(1);
                        oldNewCells.put(oldEdge, newEdge);

                        incidentCellsE = (HashSet<Cell>[]) new HashSet[2];
                        incidentCellsE[0] = new HashSet<Cell>();
                        incidentCellsE[1] = new HashSet<Cell>();

                        incidentCellsE[1].add(newFacet);

                        newIGraph.put(newEdge, incidentCellsE);

                        Cell oldEVC1 = oldEdge.getDart().getVertex();
                        Cell oldEVC2 = oldEdge.getDart().getInvolution(0).getVertex();

                        Cell newEVC1 = oldNewCells.get(oldEVC1);
                        Set<Cell>[] incidentCellsV1 = null;;

                        // System.out.println("new Edge: "+newEdge);

                        if (newEVC1 == null)
                        {
                           newEVC1 = new Cell(oldEVC1.getVertex());
                           oldNewCells.put(oldEVC1, newEVC1);

                           incidentCellsV1 = (HashSet<Cell>[]) new HashSet[2];
                           incidentCellsV1[0] = new HashSet<Cell>();
                           incidentCellsV1[1] = new HashSet<Cell>();
                           // System.out.println(newEVC1+"(vorher): "+incidentCellsV1[1]);
                           incidentCellsV1[1].add(newEdge);

                           newIGraph.put(newEVC1, incidentCellsV1);
                        }
                        else
                        {
                           incidentCellsV1 = newIGraph.get(newEVC1);
                           if (incidentCellsV1 == null)
                           {
                              incidentCellsV1 = (HashSet<Cell>[]) new HashSet[2];
                              incidentCellsV1[0] = new HashSet<Cell>();
                              incidentCellsV1[1] = new HashSet<Cell>();
                              newIGraph.put(newEVC1, incidentCellsV1);
                           }
                           // System.out.println(newEVC1+"(vorher): "+incidentCellsV1[1]);
                           incidentCellsV1[1].add(newEdge);
                        }

                        Cell newEVC2 = oldNewCells.get(oldEVC2);
                        Set<Cell>[] incidentCellsV2 = null;

                        if (newEVC2 == null)
                        {
                           newEVC2 = new Cell(oldEVC2.getVertex());
                           oldNewCells.put(oldEVC2, newEVC2);

                           incidentCellsV2 = (HashSet<Cell>[]) new HashSet[2];
                           incidentCellsV2[0] = new HashSet<Cell>();
                           incidentCellsV2[1] = new HashSet<Cell>();
                           // System.out.println(newEVC1+"(vorher): "+incidentCellsV2[1]);
                           incidentCellsV2[1].add(newEdge);

                           newIGraph.put(newEVC2, incidentCellsV2);
                        }
                        else
                        {
                           incidentCellsV2 = newIGraph.get(newEVC2);
                           if (incidentCellsV2 == null)
                           {
                              incidentCellsV2 = (HashSet<Cell>[]) new HashSet[2];
                              incidentCellsV2[0] = new HashSet<Cell>();
                              incidentCellsV2[1] = new HashSet<Cell>();
                              newIGraph.put(newEVC2, incidentCellsV2);
                           }
                           // System.out.println(newEVC2+"(vorher): "+incidentCellsV2[1]);
                           incidentCellsV2[1].add(newEdge);
                        }

                        // System.out.println(newEVC1+": "+incidentCellsV1[1]);
                        // System.out.println(newEVC2+": "+incidentCellsV2[1]);

                        incidentCellsE[0].add(newEVC1);
                        incidentCellsE[0].add(newEVC2);
                     }
                     else
                     {
                        incidentCellsE = newIGraph.get(newEdge);
                        incidentCellsE[1].add(newFacet);
                     }

                     incidentCellsF[0].add(newEdge);
                  }
                  // else
                  // {
                  // System.out.println(">>>degenerierte Kante");
                  // }
               }
            }

            // Set<Cell> cells = newIGraph.keySet();
            // Iterator<Cell> cIt = cells.iterator();
            // while (cIt.hasNext())
            // {
            // Cell curCell = cIt.next();
            // if (curCell.getDimension() == 0)
            // {
            // System.out.println(curCell.print());
            // Set<Cell>[] incidentC = newIGraph.get(curCell);
            //
            // System.out.println("\tChilds("+incidentC[0].size()+"):\n\t\t" + incidentC[0]);
            // System.out.println("\tParents("+incidentC[1].size()+"):\n\t\t" + incidentC[1]);
            // }
            // }

            System.out.println("<<<<splitted cells");
            // this seems to be correct
            Iterator<Cell> itSC = splittedCells.iterator();
            while (itSC.hasNext())
            {

               Cell oldFacet = itSC.next();
               // System.out.println(">>"+oldFacet);
               Cell newFacet = new Cell(2);
               newFacet.setData(oldFacet.getData());
               oldNewCells.put(oldFacet, newFacet);

               Set<Cell>[] incidentCellsF = (HashSet<Cell>[]) new HashSet[2];
               incidentCellsF[0] = new HashSet<Cell>();
               incidentCellsF[1] = new HashSet<Cell>();

               newIGraph.put(newFacet, incidentCellsF);
               ArrayList<Cell[]> segments = postSplitList.get(oldFacet);
               for (int i = 0; i < segments.size(); i++)
               {

                  Cell edgeV1 = segments.get(i)[0];
                  Set<Cell>[] incidentV1 = (HashSet<Cell>[]) new HashSet[2];
                  Cell edgeV2 = segments.get(i)[1];

                  Set<Cell>[] incidentV2 = (HashSet<Cell>[]) new HashSet[2];

                  if (edgeV1 != edgeV2)
                  {
                     Cell newEdge = null;
                     Set<Cell>[] incidentCellsE;
                     boolean potEdge = true;

                     if (!newIGraph.containsKey(edgeV1))
                     {
                        incidentV1[0] = new HashSet<Cell>();
                        incidentV1[1] = new HashSet<Cell>();
                        newIGraph.put(edgeV1, incidentV1);
                        potEdge = false;
                     }
                     else
                        incidentV1 = newIGraph.get(edgeV1);

                     if (!newIGraph.containsKey(edgeV2))
                     {
                        incidentV2[0] = new HashSet<Cell>();
                        incidentV2[1] = new HashSet<Cell>();
                        newIGraph.put(edgeV2, incidentV2);
                        potEdge = false;
                     }
                     else
                        incidentV2 = newIGraph.get(edgeV2);

                     if (potEdge)
                     {
                        Set<Cell> edgesAtV1 = incidentV1[1];
                        Set<Cell> edgesAtV2 = incidentV2[1];

                        Iterator<Cell> it1 = edgesAtV1.iterator();
                        while (newEdge == null && it1.hasNext())
                        {
                           Cell e1 = it1.next();
                           Iterator<Cell> it2 = edgesAtV2.iterator();
                           while (newEdge == null && it2.hasNext())
                           {
                              Cell e2 = it2.next();
                              // HOTFIX: this should not be necessary, if vertex parents are correct
                              // if (e1 == e2 && newIGraph.get(e1)[0].contains(edgeV1) &&
                              // newIGraph.get(e1)[0].contains(edgeV2))
                              // {
                              // newEdge = e1;
                              // }
                              // TODO: remove hotfix
                              if (e1 == e2)
                              {
                                 newEdge = e1;
                              }
                           }
                        }

                        // if (newEdge != null)
                        // {
                        // if (edgesAtV1.size() > 0)
                        // {
                        // Iterator<Cell> itTemp1 = edgesAtV1.iterator();
                        // while (itTemp1.hasNext())
                        // {
                        // Cell e = itTemp1.next();
                        // System.out.println("<<<<1" + e + ": " + newIGraph.get(e)[0]);
                        // }
                        // }
                        // if (edgesAtV2.size() > 0)
                        // {
                        // Iterator<Cell> itTemp2 = edgesAtV2.iterator();
                        // while (itTemp2.hasNext())
                        // {
                        // Cell e = itTemp2.next();
                        // System.out.println("<<<<2" + e + ": " + newIGraph.get(e)[0]);
                        // }
                        // }
                        // }
                     }

                     if (newEdge == null)
                     {
                        newEdge = new Cell(1);

                        incidentCellsE = (HashSet<Cell>[]) new HashSet[2];
                        incidentCellsE[0] = new HashSet<Cell>();
                        incidentCellsE[1] = new HashSet<Cell>();

                        incidentCellsE[1].add(newFacet);

                        newIGraph.put(newEdge, incidentCellsE);

                        incidentCellsE[0].add(edgeV1);
                        incidentCellsE[0].add(edgeV2);

                        incidentV1[1].add(newEdge);
                        incidentV2[1].add(newEdge);
                     }
                     else
                     {
                        incidentCellsE = newIGraph.get(newEdge);
                     }
                     incidentCellsF[0].add(newEdge);
                     incidentCellsE[1].add(newFacet);

                  }

               }
            }

            // Set<Cell> cells = newIGraph.keySet();
            // Iterator<Cell> cIt = cells.iterator();
            // while (cIt.hasNext())
            // {
            // Cell curCell = cIt.next();
            // if (curCell == null)
            // System.out.println("NUll cell in IG!!!");
            //
            // System.out.println(curCell.print());
            // Set<Cell>[] incidentC = newIGraph.get(curCell);
            //
            // System.out.println("\tChilds:\n\t\t" + incidentC[0]);
            // System.out.println("\tParents:\n\t\t" + incidentC[1]);
            // }
            //
            System.out.println("<<< incidence graph finished");
            gmo_return.buildGMapByIncidenceGraphGeneralized(newIGraph);
            System.out.println("<<< GMap finished");
         }

         return gmo_return;
      }
      else
         if (maxDimension == 3 && ncell.getDimension() == 3)
         {
            // 3D case intersection by a polyhedron
            if (vertexArrCell[0].length == 3)
            {
               for (int i = 0; i < cellVertices.size(); i++)
               {
                  double[] cellVrtxCoord = vertexArrCell[i];
                  int id = nn_Map.getIndex(cellVrtxCoord);

                  CellList incidentNCells = mapVertices.get(id).getNeighboringCells(3);
                  int j = 0;
                  while (j < incidentNCells.size() && !incidentNCells.get(j).containsPoint(cellVrtxCoord))
                  {
                     j++;
                  }

                  if (j < incidentNCells.size())
                  {
                     connectivityFunction.put(cellVertices.get(i), incidentNCells.get(j));
                  }
                  else
                  {
                     connectivityFunction.put(cellVertices.get(i), mapVertices.get(id));
                  }
               }

               for (int i = 0; i < mapVertices.size(); i++)
               {
                  double[] cellVrtxCoord = vertexArrMap[i];
                  Cell vc = mapVertices.containsPoint(vertexArrMap[i]);
                  if (ncell.containsPoint(cellVrtxCoord))
                  {
                     if (vc != null)
                     {
                        connectivityFunction.put(vc, ncell);
                     }

                  }
                  else
                  {
                     int id = nn_Cell.getIndex(cellVrtxCoord);
                     connectivityFunction.put(vc, cellVertices.get(id));
                  }
               }
               GMap gmo_return = new GMap(3);

               CellList facets = getFacetCells();
               for (int f = 0; f < facets.size(); f++)
               {

               }
               System.out.println("GMap.intersect(Cell): not implemented for dimension != 2!");
               return gmo_return;
            }
            else
            {
               System.out.println("GMap.intersect(Cell): not implemented for dimension != 2!");
               return null;
            }
         }
         else
         {
            System.out.println("GMap.intersect(Cell): not implemented for dimension != 2!");
            return null;
         }
      // TODO: 3D intersection

   }
   // private
   public void buildGMapByIncidenceGraphGeneralized(Map<Cell, Set<Cell>[]> iGraph)
   {

      Set<Cell> newCells = iGraph.keySet();
      Iterator<Cell> it = newCells.iterator();

      ArrayList<Cell> nCells = new ArrayList<Cell>();
      while (it.hasNext())
      {
         Cell cell = it.next();
         cells[cell.getDimension()].add(cell);

         if (cell.getDimension() == maxDimension)
            nCells.add(cell);
      }

      Map<Cell, ArrayList<Dart>> cellDarts = new IdentityHashMap<Cell, ArrayList<Dart>>();

      for (int c = 0; c < nCells.size(); c += 1)
      {
         this.addNCell(nCells.get(c), iGraph, cellDarts);
      }
   }

   private ArrayList<Dart> addNCell(Cell cell, Map<Cell, Set<Cell>[]> iGraph, Map<Cell, ArrayList<Dart>> allCellDarts)
   {
      if (this.maxDimension == cell.getDimension())
      {
         ArrayList<Dart> cellDarts = createCell(cell, iGraph, this.maxDimension, allCellDarts);

         // TODO:potentielle Nachbarn finden und vernähen

         for (int i = 0; i < cellDarts.size(); i++)
         {
            Dart d = cellDarts.get(i);
            // for(int j = 0; j <=d.getDimension();j++)
            int j = this.maxDimension;
            {
               if (d == d.getInvolution(j))
               {
                  ArrayList<Dart> partner = allCellDarts.get(d.getIncidentCell(j - 1));

                  for (int k = 0; k < partner.size(); k++)
                  {
                     Dart d2 = partner.get(k);
                     if (d != d2 && d2 == d2.getInvolution(j) && d.isInvolution(d2, j))
                     {
                        d.sew(d2, j);
                        break;
                     }
                  }
               }
            }
         }

         darts.addAll(cellDarts);
         return cellDarts;
      }
      return null;
   }

   private ArrayList<Dart> createCell(Cell cell, Map<Cell, Set<Cell>[]> iGraph, int maxDim,
         Map<Cell, ArrayList<Dart>> allCellDarts)
   {

      ArrayList<Dart> newDarts = new ArrayList<Dart>();

      if (cell.getDimension() == 0)
      {
         System.out.println("This should never happen ... ");
         // maybe not needed, isolated Dart

         // this.embeddings[0].add(new Cell_Embedding(embedding,true));
         Dart d = new Dart(maxDim);

         d.getIncidentCells()[0] = cell;// hier mï¿½sste ne shallowCopy rein
         // d.getIncident_cells()[0] = new Cell_Embedding(embedding,true);

         d.isCellKey(0);
         newDarts.add(d);
         if (cell.getDart() == null)
            cell.setDart(d);

      }
      else
         if (cell.getDimension() == 1)
         {
            Cell[] edgeVertexCells = new Cell[iGraph.get(cell)[0].size()];

            iGraph.get(cell)[0].toArray(edgeVertexCells);

            Dart d1 = new Dart(maxDim);
            d1.getIncidentCells()[0] = edgeVertexCells[0];// hier mï¿½sste ne shallowCopy rein
            d1.getIncidentCells()[1] = cell;// hier mï¿½sste ne shallowCopy rein
            d1.isCellKey(1);
            d1.isCellKey(1);
            newDarts.add(d1);

            if (edgeVertexCells[0].getDart() == null)
               edgeVertexCells[0].setDart(d1);

            Dart d2 = new Dart(maxDim);
            d2.getIncidentCells()[0] = edgeVertexCells[1];// hier mï¿½sste ne shallowCopy rein
            d2.getIncidentCells()[1] = cell;// hier mï¿½sste ne shallowCopy rein
            d2.isCellKey(1);
            d2.isCellKey(1);
            newDarts.add(d2);
            d1.sew(d2, 0);

            if (edgeVertexCells[1].getDart() == null)
               edgeVertexCells[1].setDart(d2);

            if (cell.getDart() == null)
               cell.setDart(d1);
         }
         else
         {
            Cell[] subCells = new Cell[iGraph.get(cell)[0].size()];

            // System.out.println(iGraph.get(cell)[0]);

            iGraph.get(cell)[0].toArray(subCells);

            // System.out.println(Arrays.toString(subCells));

            for (int sc = 0; sc < subCells.length; sc += 1)
            {
               ArrayList<Dart> recDarts = createCell(subCells[sc], iGraph, maxDim, allCellDarts);
               for (int d = 0; d < recDarts.size(); d += 1)
               {
                  recDarts.get(d).getIncidentCells()[cell.getDimension()] = cell;
                  recDarts.get(d).isCellKey(cell.getDimension());

                  for (int nd = 0; nd < newDarts.size(); nd += 1)
                  {
                     if (recDarts.get(d).isInvolution(newDarts.get(nd), cell.getDimension() - 1))
                     {
                        recDarts.get(d).sew(newDarts.get(nd), cell.getDimension() - 1);
                        break;
                     }
                  }
               }

               newDarts.addAll(recDarts);
            }
            if (cell.getDart() == null)
               cell.setDart(newDarts.get(0));
         }

      if (allCellDarts.containsKey(cell))
      {
         allCellDarts.get(cell).addAll(newDarts);
      }
      else
      {
         ArrayList<Dart> cdrts = new ArrayList<Dart>(newDarts);
         allCellDarts.put(cell, cdrts);
      }

      return newDarts;
   }

   @SuppressWarnings("unchecked")
   protected void initCellList()
   {
      cells = (ArrayList<Cell>[]) new ArrayList[maxDimension + 1];
      for (int i = 0; i <= maxDimension; i++)
      {
         cells[i] = new ArrayList<Cell>();
      }
   }

   private void ensureTIG()
   {
      if (temporaryIncidenceGraph == null)
         temporaryIncidenceGraph = new IdentityHashMap<Cell, Set<Cell>[]>();
   }
}
