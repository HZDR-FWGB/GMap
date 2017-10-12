package de.hzdr.jgm.cgeo.gmap;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Stack;

/**
 * Class for dart objects ... base element for GMaps
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class Dart
{

   /**
    * dimension
    */
   private int       dimension;
   /**
    * involutions,used in Mallet and Levy, 1999
    */
   private Dart[]    involutions;
   /**
    * incident cells, cell tuple,used in Mallet and Levy, 1999
    */
   private Cell[]    cells;
   /**
    * used for traversion in Mallet and Levy, 1999 and Brisson, 1990
    */
   private boolean   isMarked = false;
   /**
    * used in Mallet and Levy, 1999
    */
   private boolean[] isCellKey;

   /**
    * Polarity of a dart, used to enable orientability as described in Mallet (2002). A polar dart can only be sew to a
    * dart with the opposite polarity. A neutral dart can be sew to a polar dart, but its polarity changes during the
    * sew to the matching (opposite) polarity A dart with polarity = 0 is neutral, polarity =-1 is negative and polarity
    * =1 is positive.
    */
   private byte      polarity = 0;

   /**
    * needed for external searching
    */
   private boolean   visited  = false;

   /**
    * @return the visited
    */
   protected boolean isVisited()
   {
      return visited;
   }

   /**
    * visited = true
    */
   protected void visit()
   {
      this.visited = true;
   }
   
   /**
    * visited = false
    */
   protected void unVisit()
   {
      this.visited = false;
   }
   
   public Dart(int dim)
   {
      dimension = dim;
      cells = new Cell[dimension + 1];
      involutions = new Dart[dimension + 1];
      Arrays.fill(involutions, this);
      isCellKey = new boolean[dimension + 1];
      // TODO Auto-generated constructor stub
   }

   public Dart(Cell... cell_tuple)
   {
      cells = cell_tuple;
      dimension = cells.length - 1;
      involutions = new Dart[dimension + 1];
      Arrays.fill(involutions, this);
      isCellKey = new boolean[dimension + 1];
      Arrays.fill(isCellKey, true);
   }

   /**
    * copy constructor for a shallow, unsewed copy of the given dart
    * 
    * ... equal to new Dart(dart.getIncidentCells) ...
    * 
    */
   public Dart(Dart dart)
   {
      cells = dart.cells;
      involutions = new Dart[dimension + 1];
      Arrays.fill(involutions, this);
      isCellKey = new boolean[dimension + 1];
      Arrays.fill(isCellKey, true);
   }

   /**
    * sets the polarity of a dart
    * 
    * @param p
    *           polarity
    */
   public void setPolarity(int p)
   {
      if (p < 0)
         this.polarity = -1;
      else
         if (p > 0)
            this.polarity = 1;
         else
            this.polarity = 0;

      // ensure correct polarity to involutions
      // traverse all connected darts to ensure correct polarity
      Stack<Dart> st = new Stack<Dart>();

      st.push(this);

      while (st.size() > 0)
      {
         Dart dart = st.pop();
//         System.out.println(dart.print());
         byte opp = (byte) (-1 * dart.polarity);
         for (int i = 0; i <= this.dimension; i++)
         {
            if (dart != dart.involutions[i] && dart.involutions[i].polarity != opp)
            {
               if (dart.involutions[i].polarity == 0 || dart.polarity == 0)
               {
                  dart.involutions[i].polarity = opp;
                  st.push(dart.involutions[i]);
               }
            }
         }
      }
   }

   /**
    * gets the polarity of a dart
    */
   public int getPolarity()
   {
      return this.polarity;
   }

   public int getDimension()
   {
      return dimension;
   }

   public Cell[] getIncidentCells()
   {
      return cells;
   }

   public Cell getIncidentCell(int i)
   {
      return cells[i];
   }

   /**
    * return the vertex cell of this dart
    * 
    * @return
    */
   public Cell getVertex()
   {
      return cells[0];
   }

   /**
    * return the vertex coordinates of this dart
    * 
    * @return
    */
   public double[] getVertexCoordinates()
   {
      return cells[0].getVertex();
   }

   public void setIncidentCell(Cell incident_cell)
   {
      int dim = incident_cell.getDimension();
      // if (dim >= 0 && dim <= this.dimension)
      {
         cells[dim] = incident_cell;
      }
   }

   public void setIncidentCells(Cell... incident_cells)
   {
      cells = incident_cells;
   }

   public Dart[] getInvolutions()
   {
      return involutions;
   }

   public Dart getInvolution(int i)
   {
      return this.involutions[i];
   }

   /**
    * sets dim-involution bijectively
    * 
    * this.involutions[dim] = involution AND involution.involutions[dim] = this
    * 
    * Note: better use {@link sew(Dart)}
    * 
    * @param dim
    * @param involution
    */
   @SuppressWarnings("unused")
   private void setInvolution(int dim, Dart involution)
   {
      this.involutions[dim] = involution;
      involution.involutions[dim] = this;
   }

   public boolean[] getIsCellKey()
   {
      return this.isCellKey;
   }

   public void isCellKey(int dim)
   {
      this.isCellKey[dim] = true;
   }

   @SuppressWarnings("unused")
   private boolean isMarked()
   {
      return this.isMarked;
   }

   private void mark()
   {
      this.isMarked = true;
   }

   private void unmark()
   {
      this.isMarked = false;
   }

   /**
    * Evaluates this dart to be an orphan (one or more involutions reference to this dart)
    * 
    * @return lowest dimension of self-reference-involution or -1
    */
   public int isOrphan()
   {
      boolean isO = false;
      int curInv = -1;
      while (!isO && curInv < this.dimension)
      {
         curInv += 1;
         isO = this.involutions[curInv].equals(this);
      }

      if (isO)
         return curInv;
      else
         return -1;
   }

   /**
    * Checks whether another dart d2 is an involution to this dart in dimension testDim
    * 
    * @param d2
    *           second dart
    * @param testDim
    *           dimension to test
    * @return true if for both darts their embeddings only differ in testDim
    */
   public boolean isInvolution(Dart d2, int testDim)
   {
      boolean ret = false;

      ret = (this.dimension == d2.dimension);
      if (ret)
      {
         for (int d = 0; d <= this.dimension; d += 1)
         {
            if (d == testDim)
            {
               // ret = ret && !(this.incident_cells[d].equalGraphNode(d2.incident_cells[d]));
               ret = ret && !(this.cells[d] == d2.cells[d]);
            }
            else
            {
               if (this.cells[d] == null || d2.cells[d] == null)
               {
                  //
               }
               else
               {
                  // ret = ret && (this.incident_cells[d].equalGraphNode(d2.incident_cells[d]));
                  ret = ret && (this.cells[d] == d2.cells[d]);
               }
            }

            if (!ret)
               break;
         }
      }
      return ret;
   }

   /**
    * Traverse darts through a single involution (implementation from Mallet and Levy, 1999)
    * 
    * @param involutionIndex
    *           index 0<=K<=N of the involution to be used
    * @return List of traversed Darts
    * 
    */
   public ArrayList<Dart> traverse(int involutionIndex)
   {
      ArrayList<Dart> traversedDarts = new ArrayList<Dart>();
      if (involutionIndex >= 0 && involutionIndex <= this.dimension)
      {
         Stack<Dart> st = new Stack<Dart>();

         this.mark();
         st.push(this);

         while (st.size() > 0)
         {
            Dart dart = st.pop();
            traversedDarts.add(dart);
            if (!dart.involutions[involutionIndex].isMarked)
            {
               dart.involutions[involutionIndex].mark();
               st.push(dart.involutions[involutionIndex]);
            }
         }
         Dart.unmarkDartList(traversedDarts);
         return traversedDarts;
      }
      else
         return null;
   }

   /**
    * Traverse darts through a set of involutions (implementation from Mallet and Levy, 1999)
    * 
    * @param involutionIndex
    *           indices 0<=K1<=K2...<=N of the involution to be used
    * @return List of traversed Darts
    * 
    */
   public ArrayList<Dart> traverse(int[] involutionIndices)
   {
      return this.traverse(involutionIndices, false);
   }

   /**
    * Traverse darts through a set of involutions (implementation from Mallet and Levy, 1999)
    * 
    * @param involutionIndex
    *           indices 0<=K1<=K2...<=N of the involution to be used
    * @param getOrphansOnly
    *           return only orphan darts
    * 
    * @return List of traversed Darts
    * 
    */
   public ArrayList<Dart> traverse(int[] involutionIndices, boolean getOrphansOnly)
   {
      ArrayList<Dart> traversedDarts = new ArrayList<Dart>();
      ArrayList<Dart> tmpDarts = new ArrayList<Dart>();
      Stack<Dart> st = new Stack<Dart>();

      this.mark();
      st.push(this);

      while (st.size() > 0)
      {
         Dart dart = st.pop();
         if (!getOrphansOnly || (dart.isOrphan() > -1))
         {
            traversedDarts.add(dart);
         }

         tmpDarts.add(dart);

         for (int j = 0; j < involutionIndices.length; j += 1)
         {
            int involutionIndex = involutionIndices[j];
            if (involutionIndex >= 0 && involutionIndex <= this.dimension
                  && !dart.involutions[involutionIndex].isMarked)
            {
               dart.involutions[involutionIndex].mark();
               st.push(dart.involutions[involutionIndex]);
            }
         }
      }

      Dart.unmarkDartList(tmpDarts);
      return traversedDarts;
   }

   /**
    * (implementation from Mallet and Levy, 1999)
    * 
    * @param dim
    *           dimension
    * @return dart
    * 
    */
   public Dart find_cell_key(int dim)
   {
      Dart ret = null;
      int[] involutionsToBeTraversed = new int[1];
      if (dim <= 1)
         involutionsToBeTraversed[0] = dim - 1;
      else
      {
         involutionsToBeTraversed = new int[dim];
         for (int i = 0; i < dim - 1; i += 1)
         {
            involutionsToBeTraversed[i] = i;
         }
      }

      ArrayList<Dart> traversedDarts = this.traverse(involutionsToBeTraversed);
      Iterator<Dart> dit = traversedDarts.iterator();
      while (dit.hasNext())
      {
         Dart n = dit.next();
         if (n.isCellKey[dim])
         {
            return n;
         }
      }

      return ret;
   }

   /**
    * (implementation from Mallet and Levy, 1999)
    * 
    * @param embedding
    *           realization of topological object to be attached
    * 
    */
   public void dispatch_embedding(Cell embedding)
   {
      int dim = embedding.getDimension();
      int[] involutionsToBeTraversed = new int[1];
      if (dim <= 1)
         involutionsToBeTraversed[0] = dim - 1;
      else
      {
         involutionsToBeTraversed = new int[dim];
         for (int i = 0; i < dim - 1; i += 1)
         {
            involutionsToBeTraversed[i] = i;
         }
      }
      ArrayList<Dart> traversedDarts = this.traverse(involutionsToBeTraversed);
      Iterator<Dart> dit = traversedDarts.iterator();
      while (dit.hasNext())
      {
         Dart n = dit.next();
         n.cells[dim] = embedding;
      }
   }

   /**
    * (implementation from Mallet and Levy, 1999)
    * 
    * @param d2
    *           dart to be sew to this
    * @param dim
    *           sew at dim-1
    */
   public void sew(Dart d2, int dim)
   {
      this.sew(d2, dim, false);
   }

   public void sew(Dart d2, int dim, boolean debugOut)
   {
      int[] involutionsToBeTraversed = new int[1];
      if (dim <= 1)
         involutionsToBeTraversed[0] = dim - 2;
      else
      {
         involutionsToBeTraversed = new int[dim - 1];
         for (int i = 0; i <= dim - 2; i += 1)
         {
            involutionsToBeTraversed[i] = i;
            if (debugOut)
               System.out.println("Dart.sew() - Debug: Traverse involution " + i);
         }

      }

      ArrayList<Dart> tdD1 = this.traverse(involutionsToBeTraversed);
      ArrayList<Dart> tdD2 = d2.traverse(involutionsToBeTraversed);
      if (debugOut)
         System.out.println("Dart.sew() - Debug: connect involution for  " + tdD1.size() + " Darts (Dart 1) und "
               + tdD2.size() + "Darts (Dart 2)");
      for (int i = 0; i < tdD1.size(); i += 1)
      {
         int j = i;
         {
            for (int d = 0; d < dim; d += 1)
            {
               Dart k1 = tdD1.get(i).find_cell_key(d);
               Dart k2 = tdD2.get(j).find_cell_key(d);

               if (k1 != null && k2 != null && !k1.equals(k2))
               {
                  k2.isCellKey[d] = false;
               }
            }
            tdD1.get(i).involutions[dim] = tdD2.get(j);
            tdD2.get(j).involutions[dim] = tdD1.get(i);

            // this should be the correct way:
            // tdD1.get(i).updateIsOrphan();
            // tdD2.get(j).updateIsOrphan();
         }
      }

   }

   /**
    * (implementation from Mallet and Levy, 1999)
    * 
    * @param d2
    *           dart to be unsew to this
    * @param dim
    *           unsew at dim-1
    */
   public void unsew(int dim)
   {
      int[] involutionsToBeTraversed = new int[1];
      if (dim <= 1)
         involutionsToBeTraversed[0] = dim - 1;
      else
      {
         involutionsToBeTraversed = new int[dim];
         for (int i = 0; i < dim - 1; i += 1)
         {
            involutionsToBeTraversed[i] = i;
         }
      }
      ArrayList<Dart> tdD1 = this.traverse(involutionsToBeTraversed);
      for (int i = 0; i < tdD1.size(); i += 1)
      {
         Dart d2 = tdD1.get(i).involutions[dim];
         tdD1.get(i).involutions[dim] = tdD1.get(i);
         d2.involutions[dim] = d2;

      }
      // this.orphan = true;
   }

   /**
    * creates a list of indices [0 1 ... i-1 i+1 ... N]
    * 
    * @param i
    *           index to spare
    * @return List of indices without i
    */
   public int[] getInvolutionIndicesWithoutI(int i)
   {
      if (i < 0)
         i = 0;
      if (i > this.dimension)
         i = this.dimension;

      int[] alphaWithoutI = new int[this.dimension];
      int curI = 0;
      for (int j = 0; j <= this.dimension; j += 1)
      {
         if (i != j)
         {
            alphaWithoutI[curI] = j;
            curI += 1;
         }
      }

      return alphaWithoutI;
   }

   /**
    * gets all darts for i-orbit
    * 
    * @param i
    *           orbit
    * @throws IllegalArgumentException
    *            i<0 || i>this.dimension
    * @return all darts as ArrayList
    */
   public ArrayList<Dart> getIOrbit(int i)
   {
      if (i >= 0 && i <= this.dimension)
      {
         ArrayList<Dart> cDarts = this.traverse(this.getInvolutionIndicesWithoutI(i));
         return cDarts;
      }
      else
         throw new IllegalArgumentException("Illegal i (" + i + ") for i-orbit!");
   }

   public String print()
   {
      String ret = this + "\n";
      ret = ret + "\tinvolutions:\n";
      for (int i = 0; i <= dimension; i++)
         ret = ret + "\t\t" + i + ": " + involutions[i] + "\n";
      ret = ret + "\tincident cells:\n";
      for (int i = 0; i <= dimension; i++)
         ret = ret + "\t\t" + i + ": " + cells[i] + " (" + cells[i].getDimension() + ")\n";

      return ret;
   }

   // static stuff
   /**
    * sets all isMarked flags for all listed dart objects to unmarked
    * 
    * @param dLst
    *           list of dart objects
    */
   private static void unmarkDartList(ArrayList<Dart> dLst)
   {
      Iterator<Dart> it = dLst.iterator();
      while (it.hasNext())
      {
         it.next().unmark();
      }
   }
}
