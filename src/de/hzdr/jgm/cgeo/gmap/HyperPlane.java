package de.hzdr.jgm.cgeo.gmap;

import java.util.Arrays;
import de.hzdr.jgm.cgeo.Base;

/**
 * container class representing the data needed for a n-Hyperplane 
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class HyperPlane
{
   private int dimension;
   double[] normal;
   double[] position;
   
   public HyperPlane(int dimension, double[] normal, double[] position)
   {
      super();
      this.dimension = dimension;
      this.normal = normal;
      this.position = position;
   }

   public HyperPlane(int dimension)
   {
      this.dimension = dimension;
      
      normal = new double[dimension];
      position = new double[dimension];
   }
   
   public void setNormal(double... normal)
   {
      if(normal.length!=dimension)
         throw new IllegalArgumentException("Assigned normal with wrong dimension");
      
      this.normal=normal;
   }
   
   public void setPosition(double... position)
   {
      if(position.length!=dimension)
         throw new IllegalArgumentException("Assigned position with wrong dimension");
      
      this.position=position;
   }

   /**
    * @return the dimension
    */
   public int getDimension()
   {
      return dimension;
   }

   /**
    * @return the normal
    */
   public double[] getNormal()
   {
      return normal;
   }

   /**
    * @return the position
    */
   public double[] getPosition()
   {
      return position;
   }  
   
   /**
    * Evaluates, whether a is ON the hyperplain.
    *
    * @param point
    *           n-dim point coordinates
    * @return {@code TRUE} if {@code point} lies on the plain, otherwise {@code FALSE}
    */
   public boolean isOnPLane(double... point)
   {
      double[] diff = new double[dimension];
      
      for(int i = 0; i < dimension; i++)
         diff[i] = position[i]-point[i];
      
      return (Base.dot(diff, normal) == 0);
   }   
   
   public String print()
   {
      String ret = this+":\n\tnormal: "+Arrays.toString(normal)+":\n\tposition: "+Arrays.toString(position);
      return ret;
   }
}
