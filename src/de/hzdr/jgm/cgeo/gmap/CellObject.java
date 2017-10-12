package de.hzdr.jgm.cgeo.gmap;

import java.util.List;

/**
 * interface providing functionality for Celluarlar Objects
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public interface CellObject
{
   /**
    * gets the dimension
    * @return dimension
    */
   public int getDimension();
   /**
    * gets the dimension dependent measure
    * @return
    */
   public double getMeasure();
   
   /**
    * computes the circumscribed sphere
    * @return
    */
   public Sphere getCircumscribingCircle() throws NotImplementedException;
   
   /**
    * computes the center of mass by the arithmetic mean of all vertices
    * @return
    */
   public double[] getCenterOfMass();
   /**
    * gets all bordering cells
    * @return
    * @throws NotImplementedException 
    */   
   public CellList getBorderingCells() throws NotImplementedException;
   /**
    * gets all bordering cells of a particular dimension dim <this.dimension
    * @return
    * @throws NotImplementedException 
    */   
   public CellList getBorderingCells(int dim) throws NotImplementedException;
   /**
    * gets all neighboring cells
    * @return
    * @throws NotImplementedException 
    */
   public CellList getNeighboringCells() throws NotImplementedException;
   /**
    * gets all neighboring cells, connected to this cell by cells of particular size
    * @return
    * @throws NotImplementedException 
    */
   public CellList getNeighboringCells(int dim) throws NotImplementedException;
   
   /**
    * gets all vertex cells, building this cell
    */
   public CellList getVertices();
   
   /**
    * gets all vertex cells, building this cell
    */
   public List<double[]> getVertexCoordinates();   
   
}
