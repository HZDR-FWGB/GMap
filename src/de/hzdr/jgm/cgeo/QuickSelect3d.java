package de.hzdr.jgm.cgeo;

import java.util.Arrays;

/**
 * QuickSelect3d, get index (w.r.t. original array) of k-th largest element in an array of 3d double double points
 * compare wrt. height in kd-tree
 * let initial array unchanged
 *     
 *                                
 * @author teichm73
 *
 */

public class QuickSelect3d
{
   /**
    * find original index of k-th largest element wrt to height
    * arr contains actual order of elements in the array 
    * sort arr until k, take from elements original
    * investigate array only from index left until index right  
    *  
    * @param original array
    *           double[][]
    * @param arr array of indices from original
    *           int[]
    * @param k index 
    *           int          
    * @param height index in tree
    *           int           
    * @param left index
    *           int 
    * @param right index
    *           int            
    * @return int index
    */
   public static int selectKth(double[][] original, int[] arr, int k, int height, int left, int right) 
   {
      if (arr == null || arr.length <= k)
      {
         throw new Error();
      }
      int from = left;
      int to = right;

      // if from == to we reached the kth element
      while (from < to) 
      {
         int r = from; 
         int w = to;
         int mid = arr[(r + w) / 2];

         // stop if the reader and writer meets
         while (r < w) 
         {
            if (coordinatevalue(original[arr[r]], height) >= coordinatevalue(original[mid], height)) 
            { // put the large values at the end
               int tmp = arr[w];
               arr[w] = arr[r];
               arr[r] = tmp;
               w--;
            } 
            else 
            { // the value is smaller than the pivot, skip
               r++;
            }
         }

         // if we stepped up (r++) we need to step one down
         if (coordinatevalue(original[arr[r]], height) > coordinatevalue(original[mid], height))
         {
            r--;
         }
         // the r pointer is on the end of the first k elements
         if (k <= r) 
         {
            to = r;
         } 
         else 
         {
            from = r + 1;
         }
      }

      return arr[k];
   }
   
   /**
    * get value at height  
    * 
    * @param point
    *           double[]          
    * @param height index
    *           int             
    * @return double
    */
   public static double coordinatevalue(double[] point, int height)
   {
      return point[height%3];
   }
   
   public static void main(String[] args)
   {
      int n = 10; 
      double[][] original = new double[n][3]; 
      int[] arr = new int[n];
      for (int i = 0; i < n; i++) {
         arr[i] = i;
         original[i][0] = Math.random();
         original[i][1] = Math.random();
         original[i][2] = Math.random();

         System.out.println(Arrays.toString(original[i]));
      }
      
      int node = selectKth(original, arr, (int)(n/2), 0, 0, n - 1);
      System.out.println(node);
   }
}