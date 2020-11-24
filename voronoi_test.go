// MIT License: See https://github.com/pzsz/voronoi/LICENSE.md

// Author: Przemyslaw Szczepaniak (przeszczep@gmail.com)
// Port of Raymond Hill's (rhill@raymondhill.net) JavaScript implementation
// of Steven Forune's algorithm to compute Voronoi diagrams

package voronoi_test

import (
	"fmt"
	. "github.com/zzwx/voronoi"
	"math/rand"
	"testing"
)

func verifyDiagram(diagram *Diagram, edgesCount, cellsCount, perCellCount int, t *testing.T) {
	if len(diagram.Edges) != edgesCount {
		t.Errorf("Expected %d edges not %d", edgesCount, len(diagram.Edges))
	}

	if len(diagram.Cells) != cellsCount {
		t.Errorf("Expected %d cells not %d", cellsCount, len(diagram.Cells))
	}

	if perCellCount > 0 {
		for _, cell := range diagram.Cells {
			if len(cell.Halfedges) != perCellCount {
				t.Errorf("Expected per cell edge count expected %d, not %d", perCellCount, len(cell.Halfedges))
			}
		}
	}
}

func TestVoronoi2Points(t *testing.T) {
	sites := []Vertex{
		{4, 5},
		{6, 5},
	}

	diagram := ComputeDiagram(sites, NewBBox(0, 0, 10, 10), true)
	verifyDiagram(diagram, 7, 2, 4, t)
	diagram = ComputeDiagram(sites, NewBBox(0, 0, 10, 10), false)
	verifyDiagram(diagram, 1, 2, 1, t)
}

func TestVoronoi3Points(t *testing.T) {
	sites := []Vertex{
		{4, 5},
		{6, 5},
		{5, 8},
	}

	verifyDiagram(ComputeDiagram(sites, NewBBox(0, 0, 10, 10), true),
		10, 3, -1, t)
	verifyDiagram(ComputeDiagram(sites, NewBBox(0, 0, 10, 10), false),
		3, 3, 2, t)
}

func Benchmark1000(b *testing.B) {
	rand.Seed(1234567)
	b.StopTimer()
	sites := make([]Vertex, 100)
	for j := 0; j < 100; j++ {
		sites[j].X = rand.Float64() * 100
		sites[j].Y = rand.Float64() * 100
	}
	b.StartTimer()
	ComputeDiagram(sites, NewBBox(0, 0, 100, 100), true)
}

// https://github.com/jfsmig/voronoi-0/commit/ccdaaecbc3c2233212bd1f2178ca828c5c3ab440
func TestHorizontal(t *testing.T) {
	sites := make([]Vertex, 0)
	for i := 0; i < 100; i++ {
		sites = append(sites, Vertex{X: float64(i), Y: 1})
	}
	verifyDiagram(ComputeDiagram(sites, NewBBox(0, 0, 100, 100), true),
		301, 100, 4, t)
}

// https://github.com/jfsmig/voronoi-0/commit/ccdaaecbc3c2233212bd1f2178ca828c5c3ab440
func TestVertical(t *testing.T) {
	sites := make([]Vertex, 0)
	for i := 0; i < 100; i++ {
		sites = append(sites, Vertex{X: 1, Y: float64(i)})
	}
	verifyDiagram(ComputeDiagram(sites, NewBBox(0, 0, 100, 100), true),
		301, 100, 4, t)
}

// https://github.com/jfsmig/voronoi-0/commit/ccdaaecbc3c2233212bd1f2178ca828c5c3ab440
func TestSquare(t *testing.T) {
	sites := make([]Vertex, 0)
	for i := 0; i < 10; i++ {
		for j := 0; j < 10; j++ {
			sites = append(sites, Vertex{X: float64(i), Y: float64(j)})
		}
	}
	// TODO(zzwx): Confirm 0 edges correctness
	verifyDiagram(
		ComputeDiagram(sites, NewBBox(0, 10, 0, 10), true),
		0, 100, 0, t)
}

func TestSquareMore(t *testing.T) {
	sites := make([]Vertex, 0)
	for i := 0; i < 100; i++ {
		for j := 0; j < 100; j++ {
			sites = append(sites, Vertex{X: float64(i), Y: float64(j)})
		}
	}
	verifyDiagram(
		ComputeDiagram(sites, NewBBox(0, 0, 10, 10), true),
		264, 10000, 0, t)
}

func ExampleComputeDiagram_square() {
	sites := []Vertex{
		{86, 59},
		{646, 347},
		{646, 59},
		{86, 347},
	}
	bbox := NewBBox(0, 0, 800, 450)
	d := ComputeDiagram(sites, bbox, true)
	for i, cell := range d.Cells {
		var poly []Vertex
		fmt.Printf("%d:\n", i)
		if len(cell.Halfedges) > 0 {
			poly = append(poly, Vertex{X: cell.Halfedges[0].GetStartpoint().X, Y: cell.Halfedges[0].GetStartpoint().Y})
		}
		for j := 0; j < len(cell.Halfedges)-1; j++ {
			poly = append(poly, Vertex{X: cell.Halfedges[j].GetEndpoint().X, Y: cell.Halfedges[j].GetEndpoint().Y})
		}
		fmt.Printf("%+v\n", poly)
	}
	// Output:
	// 0:
	// [{X:0 Y:203} {X:366 Y:203} {X:366 Y:0} {X:0 Y:0}]
	// 1:
	// [{X:366 Y:0} {X:366 Y:203} {X:800 Y:203} {X:800 Y:0}]
	// 2:
	// [{X:366 Y:450} {X:366 Y:203} {X:0 Y:203} {X:0 Y:450}]
	// 3:
	// [{X:366 Y:203} {X:366 Y:450} {X:800 Y:450} {X:800 Y:203}]
}
