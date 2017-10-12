package de.hzdr.jgm.cgeo.gmap;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * List of cells of a particular dimension
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class CellList implements CellObject
{

   /**
    * @return
    * @see java.util.List#iterator()
    */
   public Iterator<Cell> iterator()
   {
      return cellList.iterator();
   }

   /**
    * @param o
    * @return
    * @see java.util.List#indexOf(java.lang.Object)
    */
   public int indexOf(Object o)
   {
      return cellList.indexOf(o);
   }

   private int        dimension = -1;
   private List<Cell> cellList  = new ArrayList<Cell>();

   public CellList()
   {
      // TODO Auto-generated constructor stub
   }

   public CellList(int dim)
   {
      dimension = dim;
   }

   public CellList(Collection<? extends Cell> c)
   {      
      cellList= new ArrayList<Cell>(c);
      dimension = cellList.get(0).getDimension();
   }

   public List<Cell> getList()
   {
      return cellList;
   }
   
   // delegate methods
   /**
    * @return
    * @see java.util.List#size()
    */
   public int size()
   {
      return cellList.size();
   }

   /**
    * @return
    * @see java.util.List#isEmpty()
    */
   public boolean isEmpty()
   {
      return cellList.isEmpty();
   }

   /**
    * @param o
    * @return
    * @see java.util.List#contains(java.lang.Object)
    */
   public boolean contains(Object o)
   {
      return cellList.contains(o);
   }

   /**
    * @return
    * @see java.util.List#toArray()
    */
   public Object[] toArray()
   {
      return cellList.toArray();
   }

   /**
    * @param a
    * @return
    * @see java.util.List#toArray(java.lang.Object[])
    */
   public <T> T[] toArray(T[] a)
   {
      return cellList.toArray(a);
   }

   /**
    * @param e
    * @return
    * @see java.util.List#add(java.lang.Object)
    */
   public boolean add(Cell e)
   {
      return cellList.add(e);
   }

   /**
    * @param o
    * @return
    * @see java.util.List#remove(java.lang.Object)
    */
   public boolean remove(Object o)
   {
      return cellList.remove(o);
   }

   /**
    * @param c
    * @return
    * @see java.util.List#addAll(java.util.Collection)
    */
   public boolean addAll(Collection<? extends Cell> c)
   {
      return cellList.addAll(c);
   }

   /**
    * @param index
    * @param c
    * @return
    * @see java.util.List#addAll(int, java.util.Collection)
    */
   public boolean addAll(int index, Collection<? extends Cell> c)
   {
      return cellList.addAll(index, c);
   }

   public boolean addAll(CellList otherCellList)
   {
      return cellList.addAll(otherCellList.cellList);
   }
   
   /**
    * @param c
    * @return
    * @see java.util.List#removeAll(java.util.Collection)
    */
   public boolean removeAll(Collection<?> c)
   {
      return cellList.removeAll(c);
   }

   /**
    * 
    * @see java.util.List#clear()
    */
   public void clear()
   {
      cellList.clear();
   }

   /**
    * @param index
    * @return
    * @see java.util.List#get(int)
    */
   public Cell get(int index)
   {
      return cellList.get(index);
   }

   /**
    * @param index
    * @param element
    * @see java.util.List#add(int, java.lang.Object)
    */
   public void add(int index, Cell element)
   {
      cellList.add(index, element);
   }

   /**
    * @param index
    * @return
    * @see java.util.List#remove(int)
    */
   public Cell remove(int index)
   {
      return cellList.remove(index);
   }

   // interface functionality, not overridden
   /**
    * returns the cell, that contains theses coordinates
    * 
    * @param pointCoordinates
    * @returnCell or null
    */
   public Cell containsPoint(double... pointCoordinates)
   {
      for (int i = 0; i < cellList.size(); i++)
      {
         if (cellList.get(i).containsPoint(pointCoordinates))
         {
            return cellList.get(i);
         }
      }

      return null;
   }

   public void setDimension(int dim)
   {
      dimension = dim;
   }

   
   // overrided methods
   @Override
   public int getDimension()
   {
      return dimension;
   }

   @Override
   public double getMeasure()
   {
      // TODO Auto-generated method stub
      if (dimension > 0 && dimension < 4)
      {
         double meas = 0;
         for (int c = 0; c < cellList.size(); c++)
         {
            meas = meas + cellList.get(c).getMeasure();
         }
         return meas;
      }
      return 0;
   }

   @Override
   public CellList getBorderingCells() throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public CellList getNeighboringCells() throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public CellList getBorderingCells(int dim) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public CellList getNeighboringCells(int dim) throws NotImplementedException
   {
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public CellList getVertices()
   {
      // TODO Auto-generated method stub
      CellList vCells = new CellList(0);
      for (int c = 0; c < cellList.size(); c++)
      {
         vCells.addAll(cellList.get(c).getVertices());
      }     
      return vCells;
   }

   @Override
   public List<double[]> getVertexCoordinates()
   {
      List<double[]> vrtxCoords = new ArrayList<double[]>();
      for (int c = 0; c < cellList.size(); c++)
      {
         vrtxCoords.addAll(cellList.get(c).getVertexCoordinates());
      }     
      return vrtxCoords;      
   }

   @Override
   public Sphere getCircumscribingCircle() throws NotImplementedException
   {
      // TODO Auto-generated method stub
      throw new NotImplementedException("Not implemented yet!");
   }

   @Override
   public double[] getCenterOfMass()
   {
      double[] center = new double[get(0).getDart().getVertex().getVertex().length];//null;

      //center = get(0).getCenterOfMass();
      int n = center.length;
      for (int i = 0; i < size(); i++)
      {
         double[] next = get(i).getCenterOfMass();
         for (int j = 0; j < n; j++)
         {
            center[j] = center[j] + next[j];
         }         
      }

      for (int j = 0; j < n; j++)
      {
         center[j] = center[j] / size();
      }

      return center;
   }

}
