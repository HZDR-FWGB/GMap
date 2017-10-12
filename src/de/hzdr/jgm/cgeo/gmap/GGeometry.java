package de.hzdr.jgm.cgeo.gmap;

import java.util.List;

/**
 * interface providing high-level functionality for class GMap
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public interface GGeometry
{

   /**
    * searches for a vertex cell at a particular location
    * 
    * @param pointCoordinates
    * @return Cell reference or null
    */
   public Cell searchPoint(double... pointCoordinates);

   /**
    * searches for a vertex cell at a particular location within a given tolerance
    * 
    * @param tolerance
    * @param pointCoordinates
    * @return Cell reference or null
    */
   public Cell searchPoint(double tolerance, double... pointCoordinates);

   /**
    * searches for a vertex cell at a particular location or create a new vertex cell
    * 
    * @param pointCoordinates
    * @return new or existing Cell reference
    */
   public Cell getPoint(double... pointCoordinates) throws NotImplementedException;
   /**
    * searches for a vertex cell at a particular location within a given tolerance or create a new vertex cell
    * 
    * @param tolerance
    * @param pointCoordinates
    * @return new or existing Cell reference
    */
   public Cell getPoint(double tolerance, double... pointCoordinates) throws NotImplementedException;
   
   /**
    * creates a new vertex cell at a given location
    * @param pointCoordinates
    * @return new vertex cell
    */
   public Cell addPoint(double... pointCoordinates);
   
   /**
    * adds a list of coordinates as vertices
    * @param pointCoordinates
    * @return new list of vertex cells
    */
   public CellList addPoints(double[]... pointCoordinates);
   
   /**
    * adds a list of coordinates as vertices
    * @param pointCoordinates
    * @return new list of vertex cells
    */
   public CellList addPoints(List<double[]> pointCoordinates);

   /**
    * adds a list of vertices
    * @param pointCoordinates
    * @return new list of vertex cells
    */
   public void addVertexList(List<Cell> vertexCells) throws NotImplementedException;
   
   /**
    * searches for a cell at a particular location
    * 
    * @param pointCoordinates
    * @return Cell reference or null
    */
   public Cell searchCell(double... pointCoordinates) throws NotImplementedException;

   /**
    * searches for a cell at a particular location within a given tolerance
    * 
    * @param tolerance
    * @param pointCoordinates
    * @return Cell reference or null
    */
   public Cell searchCell(double tolerance, double... pointCoordinates) throws NotImplementedException;
   
   /**
    * creates a new i-cell based on existing sequence of i-1 cells
    * @param borderingCells
    * @return cell
    */
   public Cell addCell(Cell... borderingCells);
   
   /**
    * creates a new i-cell based on existing sequence of i-1 cells
    * @param borderingCells
    * @return cell
    */
   public Cell addCell(CellList borderingCells);
   
   /**
    * finalizes cell adding by building the GMap and connecting the topology
    * @return true, when building leads to a consistent GMap
    */
   public boolean connectCells();
   
   /**
    * returns all cells for a given dimension
    * @param dimension
    * @return
    */
   public CellList getAllCells(int dimension);
   
   /**
    * returns all vertex cells
    * @param dimension
    * @return
    */
   default public CellList getVertexCells()
   {
      return getAllCells(0);
   }
   
   /**
    * returns all edge cells
    * @param dimension
    * @return
    */
   default public CellList getEdgeCells()
   {
      return getAllCells(1);
   }

   /**
    * returns all facet cells
    * @param dimension
    * @return
    */
   default public CellList getFacetCells()
   {
      return getAllCells(2);
   }
   
   /**
    * returns all polyhedron cells
    * @param dimension
    * @return
    */
   default public CellList getPolyhedronCells()
   {
      return getAllCells(3);
   }
   /**
    * creates a new n-1 GMap, representing the projected tessellation on the Plane
    * 
    * @param plane
    * @return
    */
   public GMap intersect(HyperPlane plane);
   
   /**
    * cuts this GMap along a hyperplane and splits the cutted cells
    * 
    * Note: cutted cells are removed and replaced by the splitted ones
    * 
    * @param plane
    * @return List of new cells
    */
   public CellList CutByHyperPlane(HyperPlane plane) throws NotImplementedException;
   
   /**
    * intersects this with another GMap and returns the resulting GMap
    * @param gmap another GMap
    * @return intersection GMap
    */
   public GMap intersect(GMap gmap) throws NotImplementedException;
   
   /**
    * intersects this with an n-Cell of another GMap and returns the resulting GMap
    * @param ncell of another GMap obejct
    * @return intersection GMap
    */
   public GMap intersect(Cell ncell);
   
}
