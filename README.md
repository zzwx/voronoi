[![github.com/zzwx/voronoi](doc/gobadge.svg)](https://pkg.go.dev/github.com/zzwx/voronoi)

# Voronoi Diagrams in Go

An implementation of Steven J. Fortune's algorithm to efficient Voronoi diagrams computing in Go language.
It is based on a [Raymond Hill's JavaScript implementation](https://github.com/gorhill/Javascript-Voronoi),
forked from an excellent work by Przemyslaw Szczepaniak at [github.com/pzsz/voronoi](github.com/pzsz/voronoi).

## Credits

* Forked from [github.com/pzsz/voronoi](github.com/pzsz/voronoi) and updated to the latest modifications of the original library.
* Some tests taken from [https://github.com/jfsmig/voronoi-0](https://github.com/jfsmig/voronoi-0) clone. 
* "Fortune's algorithm" by Steven J. Fortune

## Keeping up-to-date

https://github.com/gorhill/Javascript-Voronoi/commit/3fe2165aed14d3424398e04bccb120f68114eab5

## Internal differences

```
edge.lSite => edge.LeftCell.Site
edge.rSite => edge.RightCell.Site
```

```
site.voronoiId implementation:

s.cells - ordered list of cells
s.cellMap[*site] = cell (which references site)
getCell(site) returns the cell from cellMap.
``` 

## Omitted:

* `voronoiId` property
* Junkyard stuff:
  * `Voronoi.prototype.createVertex` and `this.createVertex` - we simply add a new vertex
  * `Voronoi.prototype.createHalfedge`
  * `Voronoi.prototype.recycle`
* `Voronoi.prototype.Cell.prototype.getNeighborIds`
* `Voronoi.prototype.Cell.prototype.getBbox`
* `Voronoi.prototype.Cell.prototype.pointIntersection`
* `Voronoi.prototype.quantizeSites`  
* `CREDITS.md`, `CHANGELOG.md`, `LICENSE.md`, `README.md`
* Package help
* `diagram.vertices = this.vertices;` returned to resulting Diagram.

## API Changes

Make sure to follow changed Bbox definition sequence. It is now xLeft, yTop, xRight, yBottom. 

## Usage

See tests on constructing polygons from the cells by ranging through their halfedges. 

General usage:

```go
import "github.com/zzwx/voronoi"

func useVoronoi() {
	// Sites of voronoi diagram
	sites := []voronoi.Vertex{
		voronoi.Vertex{4, 5},
		voronoi.Vertex{6, 5},
		//...
	}

	// Create bounding box
	bbox := NewBBox(0, 0, 20, 10)

	// Compute diagram and close cells (add half edges from bounding box)
	diagram := NewVoronoi().Compute(sites, bbox, true)

	// Iterate over cells & their halfedges
	for _, cell := diagram.Cells {
		for _, hedge := cell.Halfedges {
			//...
		}	
	}

	// Iterate over all edges
	for _, edge := diagram.Edge {
		//...
	}
}
```