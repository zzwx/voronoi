// Copyright 2013 Przemyslaw Szczepaniak.
// MIT License: See https://github.com/gorhill/Javascript-Voronoi/LICENSE.md

// Author: Przemyslaw Szczepaniak (przeszczep@gmail.com)
// Port of Raymond Hill's (rhill@raymondhill.net) JavaScript implementation
// of Steven  Forune's algorithm to compute Voronoi diagrams
package voronoi

import (
	"fmt"
	"sort"
)

// Cell of voronoi diagram
type Cell struct {
	// Site of the cell
	Site Vertex
	// Array of halfedges sorted counterclockwise
	Halfedges []*Halfedge
	closeMe   bool
}

func newCell(site Vertex) *Cell {
	return &Cell{Site: site}
}

func (t *Cell) prepareHalfedges() int {
	iHalfedge := len(t.Halfedges) - 1
	if t.Site.X == 9 && t.Site.Y == 10 {
		//
		fmt.Printf("")
	}

	// get rid of unused halfedges
	// rhill 2011-05-27: Keep it simple, no point here in trying
	// to be fancy: dangling edges are a typically a minority.
	for ; iHalfedge >= 0; iHalfedge-- {
		edge := t.Halfedges[iHalfedge].Edge
		if edge.Vb.Vertex == NoVertex || edge.Va.Vertex == NoVertex {
			SpliceHalfedges(&t.Halfedges, iHalfedge, 1)
			//
			//halfedges[iHalfedge] = halfedges[len(halfedges)-1]
			//halfedges = halfedges[:len(halfedges)-1]
		}
	}
	// rhill 2011-05-26: I tried to use a binary search at insertion
	// time to keep the array sorted on-the-fly (in Cell.addHalfedge()).
	// There was no real benefits in doing so, performance on
	// Firefox 3.6 was improved marginally, while performance on
	// Opera 11 was penalized marginally.
	sort.Sort(halfedgesByAngle{t.Halfedges})
	return len(t.Halfedges)
}
