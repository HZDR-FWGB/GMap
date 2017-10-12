package de.hzdr.jgm.cgeo;

import de.hzdr.jgm.cgeo.NearestNeighborfast2d;
import de.hzdr.jgm.cgeo.NearestNeighborfast3d;

/**
 * high level class for nearest-neighborhood search based on kd-tree 
 * calls from NearestNeighborfast2d and NearestNeighborfast3d
 *  
 * 
 * @author teichm73
 *
 */

public class NearestNeighbor
{
   // class attributes 
   
   public int[] heap;
   
   public double[][] original;
   
   
   // class constructors 
   
   /**
    * construct class element
    * 
    * @param originalarray
    *           double[][]
    * @return NearestNeighbor class element
    */
   public NearestNeighbor(double[][] originalarray)
   { 
      if (originalarray == null)
      {
         throw new IndexOutOfBoundsException("Array is empty!");
      }
      else
      {
         this.original = originalarray;
         int[] initialheap = new int[original.length];
         if(original[0].length == 2)
         {
            NearestNeighborfast2d.binaryheap(original, initialheap);
            this.heap = initialheap;
         }
         else if(original[0].length == 3) 
         {
            NearestNeighborfast3d.binaryheap(original, initialheap);
            this.heap = initialheap;
         }
         else
         {
            throw new IllegalArgumentException("Array dimension not supported!");
         }
      }
   }
   
   // class getters 
   
   /**
    * get NN index in original array 
    * 
    * @param point p
    *           double[]
    * @return NearestNeighbor index
    */
   public int getIndex(double[] p)
   {
      return NearestNeighborfast2d.NNindex(original, p, heap);
   }
   
   /**
    * get NN element closest to p
    * 
    * @param point p
    *           double[]
    * @return NearestNeighbor element
    */
   public double[] getElement(double[] p)
   {
      return original[this.getIndex(p)];
   }
   
   /**
    * get distance to NN
    * 
    * @param point p
    *           double[]
    * @return NearestNeighbor distance
    */
   public double getDistance(double[] p)
   {
      return Math.sqrt(NearestNeighborfast2d.squaredEuclidean2D(p, this.getElement(p)));
   }   
}

