package de.hzdr.jgm.cgeo;

import java.util.Arrays;
import de.hzdr.jgm.cgeo.QuickSelect3d;


/**
 * make a binary heap as an array representing a kd-tree for nearest neighbor search
 * store each 3d point only once
 * 
 * find index of nearest neighbor in the heap
 *  
 * 
 * @author teichm73
 *
 */

public class NearestNeighborfast3d
{
   //************************************ building the heap **********************************

   /**
    * built kd-tree from double[][] original and store element index in heap 
    * represents a complete binary tree / binary heap
    * transforms 0,...,8 to 5,3,7,1,4,6,8,0,2
    * 
    * order of the heap is given by                0
    *                                        1           2
    *                                     3     4     5     6
    *                                   7   8  ... 
    * 
    * @param original array
    *           double[][]
    * @param arr array
    *           int[]
    * @param heap array, same length as arr
    *           int[]    
    * @return void
    */
   public static void binaryheap(double[][] original, int[] heap)
   {
      if (original.length > 0)
      {
         int[] arr = new int[original.length]; // initial order of element in original
         Arrays.setAll(arr, i -> i);
         binaryheapInPlace(original, arr, heap, 0, 0, 0, original.length - 1); // sort all
      }
   }
   
   /**
    * heap iterated 
    * 
    * @param original array
    *           double[][]
    * @param arr array
    *           int[]
    * @param heap array
    *           int[]                     
    * @param node index
    *           int  
    * @param height index
    *           int   
    * @param left index
    *           int 
    * @param right index
    *           int            
    * @return void
    */
   private static void binaryheapInPlace(double[][] original, int[] arr, int[] heap, int node, int height, int left, int right)
   {
      if (left > right)
      {
         return;
      }
      int index = left + levelmid(right - left + 1); // Maximum number of single-direction edges in leveled binary tree with n nodes, see http://oeis.org/A279521
      heap[node] = QuickSelect3d.selectKth(original, arr, index, height, left, right); // sort arr until k-th largest element and add to heap

      if(2*node + 1 < arr.length)
      {
         binaryheapInPlace(original, arr, heap, 2*node + 1, height + 1 , left, index - 1); // left branch
         if(2*node + 2 < arr.length)
         {
            binaryheapInPlace(original, arr, heap, 2*node + 2 , height + 1, index + 1, right); // right branch     
         }   
      }
   }
   /**
    * get index to obtain complete binary tree 
    *  
    * @param n length
    *           int          
    * @return int middle index
    */
   private static int levelmid(int n)
   {
      int maxheight = (int)(Math.log(n)/Math.log(2)); // height of the tree of n elements
      return Math.min((int)Math.pow(2, maxheight) - 1, n - (int)Math.pow(2, maxheight - 1));
   }

   //************************************ methods on the heap **********************************

   /**
    * go down the tree until leaf to find a candiate index for the nearest neighbor of point 
    * starting in node at height
    * ignore any marks
    * 
    * @param original array
    *           double[][]
    * @param point to find NN for
    *           double[]  
    * @param heap array
    *           int[]                     
    * @param node index
    *           int  
    * @param height index
    *           int            
    * @return int index in the heap
    */
   public static int NNcandidateindex(double[][] original, double[] point, int[] heap, int node, int height)
   {
      int index = 0;
      if (node < heap.length) // finished in leaf
      {
         if (QuickSelect3d.coordinatevalue(point,height) <= QuickSelect3d.coordinatevalue(original[heap[node]],height))
         {
            index = NNcandidateindex(original, point, heap, 2*node + 1, height + 1); // left branch
         }
         else
         {
            index = NNcandidateindex(original, point, heap, 2*node + 2 , height + 1); // right branch    
         }
      }
      else
      {
         index = (int)((node-1)/2); // no nodes under that
      }
      return index; // heap[index] is the candidate
   }  
   
   /**
    * Euclidean distance between two double double points 
    * 
    * @param a point 
    *           double[]  
    * @param b point
    *           double[]                     
    * @return double distance
    */
   public static double squaredEuclidean3D(double[] a, double[] b)
   {
      double ax = a[0]; // x coord
      double ay = a[1]; // y coord  
      double az = a[2]; // z coord 
      
      double bx = b[0]; // x coord
      double by = b[1]; // y coord
      double bz = b[2]; // z coord 

      double d = (ax - bx)*(ax - bx) + (ay - by)*(ay - by) + (az - bz)*(az - bz);

      return d;
   }
   
   /**
    * find the index of the overall nearest neighbor element in original of point p
    *  
    * @param original array
    *           double[][]
    * @param point to find NN for
    *           double[]  
    * @param heap array
    *           int[]                      
    * @return int node index
    */
   public static int NNindex(double[][] original, double[] point, int[] heap)
   {
      double[] NN = new double[2]; // initialize with root
      NN[0] = 0; // current best
      NN[1] = squaredEuclidean3D(point, original[heap[0]]);// initial NN distance      
      if (NN[1] == 0.0) // check for immediate stop
      {
         return heap[(int) NN[0]];
      }
      // get candidate leaf index
      int startnode = NNcandidateindex(original, point, heap, 0, 0); // initial leaf
      NNindexrecursive(original, point, heap, 0, startnode, NN); // recursion
  
      return heap[(int) NN[0]];
   }
   
   /**
    * go up the tree starting in node at height
    * recursion method
    * 
    * @param original array
    *           double[][]
    * @param point to find NN for
    *           double[]  
    * @param heap array
    *           int[]  
    * @param originalnode index, from that I stated
    *           int
    * @param node index
    *           int  
    * @param NN array of current best node and its squared distance
    *           int[]                               
    * @return void
    */
   private static void NNindexrecursive(double[][] original, double[] point, int[] heap, int originalnode, int node, double[] NN)
   {
      int height = (int)(Math.log(node + 1) / Math.log(2)); // height of node
      while (node > originalnode) // not finished in root  (do not replace > by >= !!!)
      { 
         //System.out.println("NNindex " + node + " " + heap[node]);
         double newdistance = squaredEuclidean3D(point, original[heap[node]]);
         if (newdistance < NN[1])
         {
            NN[0] = node; 
            NN[1] = newdistance;    
         }
         if (NN[1] == 0.0) // check for immediate stop
         {
            return;
         }
         // check if brother node exists and for crossing of the split plane
         if (node%2 == 1) // current node is a left child
         {
            if (node + 1 < heap.length && splitplanecrossing(original, point, heap, (int)((node - 1)/2), height - 1, NN[1])) // right child does exist and crossing occours
            {
               int startnode = NNcandidateindex(original, point, heap, node + 1, height); // new start node
               NNindexrecursive(original, point, heap, node + 1, startnode, NN); 
            }
         }
         else
         {
            if (splitplanecrossing(original, point, heap, (int)((node - 1)/2), height - 1, NN[1])) // check plane crossing
            {
               int startnode = NNcandidateindex(original, point, heap, node - 1, height); // new star node
               NNindexrecursive(original, point, heap, node - 1, startnode, NN);
            }
         }
         node = (int)((node - 1)/2); // go up in the tree
         height--;   
      }
      // check original node
      //System.out.println("NNindex " + node + " " + heap[node]);
      double newdistance = squaredEuclidean3D(point, original[heap[node]]);
      if (newdistance < NN[1])
      {
         NN[0] = node; 
         NN[1] = newdistance;    
      }
      if (NN[1] <= 0.0) // check for immediate stop
      {
         return;
      }
   } 
   
   /**
    * decide whether the splitting plane trough node and p 
    * intersects the sphere with radius NNdistance around p  
    * decides when to go down the other subtree 
    * 
    * @param original array
    *           double[][]
    * @param point to find NN for
    *           double[]  
    * @param heap array
    *           int[]  
    * @param node to be checked
    *           int 
    * @param height current height
    *           int  
    * @param NNdistance current squared NN dist 
    *           double                                         
    * @return boolean
    */
   private static boolean splitplanecrossing(double[][] original, double[] point, int[] heap, int node, int height, double NNdistance)
   {
      double sqdist = QuickSelect3d.coordinatevalue(point, height) - QuickSelect3d.coordinatevalue(original[heap[node]], height);
      boolean b = sqdist*sqdist < NNdistance;
      return b;
   }
   
   public static void main(String[] args)
   {
      int n = 1000000; 
      double[][] original = new double[n][3]; 
      for (int i = 0; i < n; i++) {
         original[i][0] = Math.random();
         original[i][1] = Math.random();
         original[i][2] = Math.random();
      }

    
      double[] p = {0.3, 0.1, 0.5};
      original[0] = p;
      int[] heap = new int[original.length]; 
      binaryheap(original, heap);
      System.out.println(Arrays.toString(original[NNindex(original, p, heap)]));
   }
}