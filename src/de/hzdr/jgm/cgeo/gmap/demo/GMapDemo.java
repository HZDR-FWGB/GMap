/**
 * 
 */
package de.hzdr.jgm.cgeo.gmap.demo;

import java.util.ArrayList;
import de.hzdr.jgm.cgeo.gmap.Cell;
import de.hzdr.jgm.cgeo.gmap.CellList;
import de.hzdr.jgm.cgeo.gmap.Dart;
import de.hzdr.jgm.cgeo.gmap.GMap;

/**
 * Simple demo for GMap package
 * 
 * @author P. Menzel - Helmholtz-Institut Freiberg for Resource Technology
 *
 */
public class GMapDemo
{
   /**
    * 
    */
   public GMapDemo()
   {
      // TODO Auto-generated constructor stub
   }

   public GMap exampleSimple3()
   {
      System.out.println("====exampleSimple3()===");
      long start = System.currentTimeMillis(); 
      
      ArrayList<double[]> points = new ArrayList<double[]>();
      points.add(new double[]{0,1});
      points.add(new double[]{0,0});
      points.add(new double[]{1,0});
      points.add(new double[]{1,1});
      points.add(new double[]{0.5,1.5});
      
      ArrayList<int[]> polys = new ArrayList<int[]>();
      polys.add(new int[]{0,1,2});
      polys.add(new int[]{0,2,3});
      polys.add(new int[]{0,3,4});
      
      GMap simple_gmo = GMap.createPolygonMesh(points, polys);
      
      long diff = System.currentTimeMillis() - start;
      System.out.println("Building as high-level poly-mesh  "
            + diff+ "ms");
      
      simple_gmo.check();
      
      CellList vertices = simple_gmo.getVertexCells();  
      System.out.println(vertices.get(0).print());
      System.out.println("Cumulated point measure ("+vertices.size()+" vertices): "+vertices.getMeasure());

      CellList edges = simple_gmo.getEdgeCells();  
      System.out.println(edges.get(0).print() + "; "+edges.get(0).getMeasure());      
      System.out.println("Cumulated edge lengths ("+edges.size()+" edges): "+edges.getMeasure());   
      
      CellList facets = simple_gmo.getFacetCells();   
      System.out.println(facets.get(0).print() + "; "+facets.get(0).getMeasure());
      System.out.println("Cumulated facet area ("+facets.size()+" facets): "+facets.getMeasure()); 
      
      return simple_gmo;
   }
   
   public void exampleSimple2()
   {
      System.out.println("====exampleSimple2()===");
      long start = System.currentTimeMillis(); 
      Cell v1 = new Cell(0,1,0);
      Cell v2 = new Cell(0,0,0);
      Cell v3 = new Cell(1,0,0);
      Cell v4 = new Cell(1,1,0);
      Cell v5 = new Cell(0.5,1.5,0);      
      
      Cell e1 = new Cell(1);
      Cell e2 = new Cell(1);
      Cell e3 = new Cell(1);
      Cell e4 = new Cell(1);
      Cell e5 = new Cell(1);
      Cell e6 = new Cell(1);
      Cell e7 = new Cell(1);
      
      Cell f1 = new Cell(2);
      Cell f2 = new Cell(2);
      Cell f3 = new Cell(2);
      
      //triangle 1
      Dart d11 = new Dart(v1,e1,f1);
      Dart d12 = new Dart(v2,e1,f1);
      Dart d13 = new Dart(v2,e2,f1);
      Dart d14 = new Dart(v3,e2,f1);
      Dart d15 = new Dart(v3,e4,f1);
      Dart d16 = new Dart(v1,e4,f1);
      
    d11.sew(d12, 0);
    d12.sew(d13, 1);
    d13.sew(d14, 0);
    d14.sew(d15, 1);
    d15.sew(d16, 0);
    d16.sew(d11, 1);
    
    Dart d21 = new Dart(v1,e4,f2);
    Dart d22 = new Dart(v3,e4,f2);
    Dart d23 = new Dart(v3,e3,f2);
    Dart d24 = new Dart(v4,e3,f2);
    Dart d25 = new Dart(v4,e5,f2);
    Dart d26 = new Dart(v1,e5,f2);
    
    d21.sew(d22, 0);
    d22.sew(d23, 1);
    d23.sew(d24, 0);
    d24.sew(d25, 1);
    d25.sew(d26, 0);
    d26.sew(d21, 1);      
    
    Dart d31 = new Dart(v1,e5,f3);
    Dart d32 = new Dart(v4,e5,f3);
    Dart d33 = new Dart(v4,e7,f3);
    Dart d34 = new Dart(v5,e7,f3);
    Dart d35 = new Dart(v5,e6,f3);
    Dart d36 = new Dart(v1,e6,f3);
    
    d31.sew(d32, 0);
    d32.sew(d33, 1);
    d33.sew(d34, 0);
    d34.sew(d35, 1);
    d35.sew(d36, 0);
    d36.sew(d31, 1);  
    
    d15.sew(d22, 2);
    d16.sew(d21, 2);
    d26.sew(d31, 2);
    d25.sew(d32, 2);      
    
    GMap simple_gmo = new GMap(2);
    simple_gmo.addDart(d11);
    simple_gmo.addDart(d12);
    simple_gmo.addDart(d13);
    simple_gmo.addDart(d14);
    simple_gmo.addDart(d15);
    simple_gmo.addDart(d16);
    simple_gmo.addDart(d21);
    simple_gmo.addDart(d22);
    simple_gmo.addDart(d23);
    simple_gmo.addDart(d24);
    simple_gmo.addDart(d25);
    simple_gmo.addDart(d26);
    simple_gmo.addDart(d31);
    simple_gmo.addDart(d32);
    simple_gmo.addDart(d33);
    simple_gmo.addDart(d34);
    simple_gmo.addDart(d35);
    simple_gmo.addDart(d36);

    long diff = System.currentTimeMillis() - start;
    System.out.println("Building with addDart without searching "
          + diff+ "ms");
    
    simple_gmo.check();
   
    CellList vertices = simple_gmo.getVertexCells();  
    System.out.println(vertices.get(0).print());
    System.out.println("Cumulated point measure ("+vertices.size()+" vertices): "+vertices.getMeasure());

    CellList edges = simple_gmo.getEdgeCells();  
    System.out.println(edges.get(0).print() + "; "+edges.get(0).getMeasure());      
    System.out.println("Cumulated edge lengths ("+edges.size()+" edges): "+edges.getMeasure());   
    
    CellList facets = simple_gmo.getFacetCells();   
    System.out.println(facets.get(0).print() + "; "+facets.get(0).getMeasure());
    System.out.println("Cumulated facet area ("+facets.size()+" facets): "+facets.getMeasure());     
    
   }
   public void exampleSimple1()
   {      
      //create cells
      System.out.println("====exampleSimple1()===");
      long start = System.currentTimeMillis(); 
      
      Cell v1 = new Cell(0,1);
      Cell v2 = new Cell(0,0);
      Cell v3 = new Cell(1,0);
      Cell v4 = new Cell(1,1);
      Cell v5 = new Cell(0.5,1.5);      
      
      Cell e1 = new Cell(1);
      Cell e2 = new Cell(1);
      Cell e3 = new Cell(1);
      Cell e4 = new Cell(1);
      Cell e5 = new Cell(1);
      Cell e6 = new Cell(1);
      Cell e7 = new Cell(1);
      
      Cell f1 = new Cell(2);
      Cell f2 = new Cell(2);
      Cell f3 = new Cell(2);
      
      //triangle 1
      Dart d11 = new Dart(v1,e1,f1);
      Dart d12 = new Dart(v2,e1,f1);
      Dart d13 = new Dart(v2,e2,f1);
      Dart d14 = new Dart(v3,e2,f1);
      Dart d15 = new Dart(v3,e4,f1);
      Dart d16 = new Dart(v1,e4,f1);
      
      Dart d21 = new Dart(v1,e4,f2);
      Dart d22 = new Dart(v3,e4,f2);
      Dart d23 = new Dart(v3,e3,f2);
      Dart d24 = new Dart(v4,e3,f2);
      Dart d25 = new Dart(v4,e5,f2);
      Dart d26 = new Dart(v1,e5,f2);    
      
      Dart d31 = new Dart(v1,e5,f3);
      Dart d32 = new Dart(v4,e5,f3);
      Dart d33 = new Dart(v4,e7,f3);
      Dart d34 = new Dart(v5,e7,f3);
      Dart d35 = new Dart(v5,e6,f3);
      Dart d36 = new Dart(v1,e6,f3);
      
      GMap simple_gmo = new GMap(2);
      simple_gmo.addDart(d11,true);
      simple_gmo.addDart(d12,true);
      simple_gmo.addDart(d13,true);
      simple_gmo.addDart(d14,true);
      simple_gmo.addDart(d15,true);
      simple_gmo.addDart(d16,true);
      simple_gmo.addDart(d21,true);
      simple_gmo.addDart(d22,true);
      simple_gmo.addDart(d23,true);
      simple_gmo.addDart(d24,true);
      simple_gmo.addDart(d25,true);
      simple_gmo.addDart(d26,true);
      simple_gmo.addDart(d31,true);
      simple_gmo.addDart(d32,true);
      simple_gmo.addDart(d33,true);
      simple_gmo.addDart(d34,true);
      simple_gmo.addDart(d35,true);
      simple_gmo.addDart(d36,true);
      
      long diff = System.currentTimeMillis() - start;
      System.out.println("Building with addDart WITH searching "
            + diff+ "ms");
      
      simple_gmo.check();
      
      CellList vertices = simple_gmo.getVertexCells();  
      System.out.println(vertices.get(0).print());
      System.out.println("Cumulated point measure ("+vertices.size()+" vertices): "+vertices.getMeasure());

      CellList edges = simple_gmo.getEdgeCells();  
      System.out.println(edges.get(0).print() + "; "+edges.get(0).getMeasure());      
      System.out.println("Cumulated edge lengths ("+edges.size()+" edges): "+edges.getMeasure());   
      
      CellList facets = simple_gmo.getFacetCells();   
      System.out.println(facets.get(0).print() + "; "+facets.get(0).getMeasure());
      System.out.println("Cumulated facet area ("+facets.size()+" facets): "+facets.getMeasure());     

    
   }

   private void exampleSimple3D2()
   {
      System.out.println("====exampleSimple3D2()===");
      
      ArrayList<double[]> points = new ArrayList<double[]>();
      points.add(new double[]{0,0,0});
      points.add(new double[]{1,0,0});
      points.add(new double[]{0.5,1,0});
      points.add(new double[]{0.5,0.5,1});
      points.add(new double[]{1.5,2,0.5});
      
      ArrayList<int[]> polys = new ArrayList<int[]>();
      polys.add(new int[]{0,1,2,3});
      polys.add(new int[]{1,3,2,4});
      
      GMap gmo = GMap.createTetrahedronMesh(points, polys);

      if(gmo.check())
      {
      CellList facets = gmo.getFacetCells();
      for(int i = 0; i < facets.size();i++)
         System.out.println("Area facet "+i+": "+facets.get(i).getMeasure());
      
      System.out.println("Volume 1st tetrahedron: "+gmo.getPolyhedronCells().get(0).getMeasure());
      System.out.println("Volume 2nd tetrahedron: "+gmo.getPolyhedronCells().get(1).getMeasure());
      }
   }
   
   private void exampleSimple3D()
   {
      System.out.println("====exampleSimple3D()===");
      
      //TODO build 2 connected tetrahedrons
      //1) 1,2,3,4
      //2) 2,3,4,5
      long start = System.currentTimeMillis(); 
      GMap gmo = new GMap(3);
      
      Cell v1 = gmo.addPoint(0, 0, 0);
      Cell v2 = gmo.addPoint(1, 0, 0);
      Cell v3 = gmo.addPoint(0.5, 1, 0);
      Cell v4 = gmo.addPoint(0.5, 0.5, 1);
      Cell v5 = gmo.addPoint(1.5, 2, 0.5);      
      
      Cell e1 = gmo.addCell(v1,v2);//new Cell(1);//v1,v2
      Cell e2 = gmo.addCell(v2,v3);//v2,v3
      Cell e3 = gmo.addCell(v3,v1);//v3,v1
      Cell e4 = gmo.addCell(v1,v4);//v1,v4
      Cell e5 = gmo.addCell(v2,v4);//v2,v4
      Cell e6 = gmo.addCell(v3,v4);//v3,v4
      Cell e7 = gmo.addCell(v2,v5);//v2,v5
      Cell e8 = gmo.addCell(v3,v5);//v3,v5
      Cell e9 = gmo.addCell(v4,v5);//v4,v5
      
      Cell f1 = gmo.addCell(e1,e2,e3);//e1,e2,e3
      Cell f2 = gmo.addCell(e1,e4,e5);//e1,e4,e5
      Cell f3 = gmo.addCell(e2,e5,e6);//e2,e5,e6
      Cell f4 = gmo.addCell(e3,e6,e4);//e3,e6,e4
      Cell f5 = gmo.addCell(e5,e7,e9);//e5,e7,e9
      Cell f6 = gmo.addCell(e2,e7,e8);//e2,e7,e8
      Cell f7 = gmo.addCell(e6,e8,e9);//e6,e8,e9
      
      Cell t1 = gmo.addCell(f1,f2,f3,f4);
      Cell t2 = gmo.addCell(f3,f5,f6,f7);
      
      gmo.connectCells();
      long diff = System.currentTimeMillis() - start;
      System.out.println("Building tetrahedral mesh by add: "
            + diff+ "ms");
      
      System.out.println("Area facet 1: "+f1.getMeasure());
      System.out.println("Area facet 2: "+f2.getMeasure());
      System.out.println("Area facet 3: "+f3.getMeasure());
      System.out.println("Area facet 4: "+f4.getMeasure());
      System.out.println("Area facet 5: "+f5.getMeasure());
      System.out.println("Area facet 6: "+f6.getMeasure());
      System.out.println("Area facet 7: "+f7.getMeasure());
      
      System.out.println("Volume 1st tetrahedron: "+t1.getMeasure());
      System.out.println("Volume 2nd tetrahedron: "+t2.getMeasure());
   }


   private void exampleIntersect2D()
   {
      System.out.println("====exampleIntersect2D()===");
      GMap rect_gmo = GMap.createQuad(new double[]{-0.5,-0.5},new double[]{.5,.5});      
      
      System.out.println(rect_gmo.getFacetCells().getMeasure());

      ArrayList<double[]> points = new ArrayList<double[]>();
      points.add(new double[]{0,0});
      points.add(new double[]{-1.5,0});
      points.add(new double[]{0,1.5});
      points.add(new double[]{1.5,0});
      points.add(new double[]{0,-1.5});
      
      ArrayList<int[]> polys = new ArrayList<int[]>();
      polys.add(new int[]{0,1,2});
      polys.add(new int[]{0,2,3});
      polys.add(new int[]{0,3,4});
      polys.add(new int[]{0,4,1});
      
      GMap gmo = GMap.createPolygonMesh(points, polys);
      
      gmo.check();
      
      System.out.println("test intersetion - simple");
      gmo.intersect(rect_gmo.getFacetCells().get(0));
     
   }

   private void exampleIntersect2Dv2()
   {
      System.out.println("====exampleIntersect2Dv2()===");
      GMap rect_gmo = GMap.createQuad(new double[]{-0.5,-0.5},new double[]{.5,.5});      
      
      System.out.println(rect_gmo.getFacetCells().getMeasure());

      ArrayList<double[]> points = new ArrayList<double[]>();
      points.add(new double[]{-0.1,0});
      points.add(new double[]{0,-0.1});
      points.add(new double[]{0.1,0});
      points.add(new double[]{0,0.1});
      
      points.add(new double[]{-1.5,0});
      points.add(new double[]{0,1.5});
      points.add(new double[]{1.5,0});
      points.add(new double[]{0,-1.5});
      
      ArrayList<int[]> polys = new ArrayList<int[]>();
      polys.add(new int[]{0,1,2,3});
      polys.add(new int[]{0,3,5,4});
      polys.add(new int[]{3,2,6,5});
      polys.add(new int[]{2,1,7,6});
      polys.add(new int[]{1,0,4,7});
      
      GMap gmo = GMap.createPolygonMesh(points, polys);
      
      gmo.check();
      
      System.out.println("test intersetion - with covered cells");           
      GMap cutGmo = gmo.intersect(rect_gmo.getFacetCells().get(0));
      cutGmo.check();
   }
   
   
   
   /**
    * @param args
    */
   public static void main(String[] args)
   {
      GMapDemo  demo = new GMapDemo();
      
      demo.exampleSimple1();
      demo.exampleSimple2();
      demo.exampleSimple3();        
      demo.exampleSimple3D();
      demo.exampleSimple3D2();      
      demo.exampleIntersect2Dv2();
//           
   }




}
