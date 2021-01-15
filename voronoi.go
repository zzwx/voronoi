// MIT License: See https://github.com/pzsz/voronoi/LICENSE.md

// Author: Przemyslaw Szczepaniak (przeszczep@gmail.com)
// Port of Raymond Hill's (rhill@raymondhill.net) JavaScript implementation
// of Steven Forune's algorithm to compute Voronoi diagrams
package voronoi

import (
	"fmt"
	"math"
	"sort"
)

type Voronoi struct {
	cells []*Cell
	edges []*Edge

	cellsMap map[Vertex]*Cell

	beachline        rbTree
	circleEvents     rbTree
	firstCircleEvent *circleEvent
}

type Diagram struct {
	Cells []*Cell
	Edges []*Edge
	//	EdgesVertices map[Vertex]EdgeVertex
}

func (s *Voronoi) getCell(site Vertex) *Cell {
	ret := s.cellsMap[site]
	if ret == nil {
		panic(fmt.Sprintf("Couldn't find cell for site %v", site))
	}
	return ret
}

//
func (s *Voronoi) createEdge(lSite, rSite *Cell, va, vb Vertex) *Edge {
	edge := newEdge(lSite, rSite)

	s.edges = append(s.edges, edge)

	if va != NoVertex {
		s.setEdgeStartpoint(edge, lSite, rSite, va)
	}
	if vb != NoVertex {
		s.setEdgeEndpoint(edge, lSite, rSite, vb)
	}
	lSite.Halfedges = append(lSite.Halfedges, newHalfedge(edge, lSite, rSite))
	rSite.Halfedges = append(rSite.Halfedges, newHalfedge(edge, rSite, lSite))
	return edge
}

func (s *Voronoi) createBorderEdge(LeftCell *Cell, va, vb Vertex) *Edge {
	edge := newEdge(LeftCell, nil)
	edge.Va.Vertex = va
	edge.Vb.Vertex = vb

	s.edges = append(s.edges, edge)
	return edge
}

func (s *Voronoi) setEdgeStartpoint(edge *Edge, LeftCell, RightCell *Cell, vertex Vertex) {
	if edge.Va.Vertex == NoVertex && edge.Vb.Vertex == NoVertex {
		edge.Va.Vertex = vertex
		edge.LeftCell = LeftCell
		edge.RightCell = RightCell
	} else if edge.LeftCell == RightCell {
		edge.Vb.Vertex = vertex
	} else {
		edge.Va.Vertex = vertex
	}
}

func (s *Voronoi) setEdgeEndpoint(edge *Edge, LeftCell, RightCell *Cell, vertex Vertex) {
	s.setEdgeStartpoint(edge, RightCell, LeftCell, vertex)
}

type Beachsection struct {
	node        *rbNode
	site        Vertex
	circleEvent *circleEvent
	edge        *Edge
}

// rbNodeValue intergface
func (s *Beachsection) bindToNode(node *rbNode) {
	s.node = node
}

// rbNodeValue intergface
func (s *Beachsection) getNode() *rbNode {
	return s.node
}

// calculate the left break point of a particular beach section,
// given a particular sweep line
func leftBreakPoint(arc *Beachsection, directrix float64) float64 {
	// http://en.wikipedia.org/wiki/Parabola
	// http://en.wikipedia.org/wiki/Quadratic_equation
	// h1 = x1,
	// k1 = (y1+directrix)/2,
	// h2 = x2,
	// k2 = (y2+directrix)/2,
	// p1 = k1-directrix,
	// a1 = 1/(4*p1),
	// b1 = -h1/(2*p1),
	// c1 = h1*h1/(4*p1)+k1,
	// p2 = k2-directrix,
	// a2 = 1/(4*p2),
	// b2 = -h2/(2*p2),
	// c2 = h2*h2/(4*p2)+k2,
	// x = (-(b2-b1) + Math.sqrt((b2-b1)*(b2-b1) - 4*(a2-a1)*(c2-c1))) / (2*(a2-a1))
	// When x1 become the x-origin:
	// h1 = 0,
	// k1 = (y1+directrix)/2,
	// h2 = x2-x1,
	// k2 = (y2+directrix)/2,
	// p1 = k1-directrix,
	// a1 = 1/(4*p1),
	// b1 = 0,
	// c1 = k1,
	// p2 = k2-directrix,
	// a2 = 1/(4*p2),
	// b2 = -h2/(2*p2),
	// c2 = h2*h2/(4*p2)+k2,
	// x = (-b2 + Math.sqrt(b2*b2 - 4*(a2-a1)*(c2-k1))) / (2*(a2-a1)) + x1

	// change code below at your own risk: care has been taken to
	// reduce errors due to computers' finite arithmetic precision.
	// Maybe can still be improved, will see if any more of this
	// kind of errors pop up again.
	site := arc.site
	rfocx := site.X
	rfocy := site.Y
	pby2 := rfocy - directrix
	// parabola in degenerate case where focus is on directrix
	if pby2 == 0 {
		return rfocx
	}

	lArc := arc.getNode().previous
	if lArc == nil {
		return math.Inf(-1)
	}
	site = lArc.value.(*Beachsection).site
	lfocx := site.X
	lfocy := site.Y
	plby2 := lfocy - directrix
	// parabola in degenerate case where focus is on directrix
	if plby2 == 0 {
		return lfocx
	}
	hl := lfocx - rfocx
	aby2 := 1/pby2 - 1/plby2
	b := hl / plby2
	if aby2 != 0 {
		return (-b+math.Sqrt(b*b-2*aby2*(hl*hl/(-2*plby2)-lfocy+plby2/2+rfocy-pby2/2)))/aby2 + rfocx
	}
	// both parabolas have same distance to directrix, thus break point is midway
	return (rfocx + lfocx) / 2
}

// calculate the right break point of a particular beach section,
// given a particular directrix
func rightBreakPoint(arc *Beachsection, directrix float64) float64 {
	rArc := arc.getNode().next
	if rArc != nil {
		return leftBreakPoint(rArc.value.(*Beachsection), directrix)
	}
	site := arc.site
	if site.Y == directrix {
		return site.X
	}
	return math.Inf(1)
}

func (s *Voronoi) detachBeachsection(beachsection *Beachsection) {
	s.detachCircleEvent(beachsection)         // detach potentially attached circle event
	s.beachline.removeNode(beachsection.node) // remove from RB-tree
}

type BeachsectionPtrs []*Beachsection

func (s *BeachsectionPtrs) appendLeft(b *Beachsection) {
	*s = append(*s, b)
	for id := len(*s) - 1; id > 0; id-- {
		(*s)[id] = (*s)[id-1]
	}
	(*s)[0] = b
}

func (s *BeachsectionPtrs) appendRight(b *Beachsection) {
	*s = append(*s, b)
}

func (s *Voronoi) removeBeachsection(beachsection *Beachsection) {
	circle := beachsection.circleEvent
	x := circle.x
	y := circle.ycenter
	vertex := Vertex{x, y}
	previous := beachsection.node.previous
	next := beachsection.node.next
	disappearingTransitions := BeachsectionPtrs{beachsection}
	abs_fn := math.Abs

	// remove collapsed beachsection from beachline
	s.detachBeachsection(beachsection)

	// there could be more than one empty arc at the deletion point, this
	// happens when more than two edges are linked by the same vertex,
	// so we will collect all those edges by looking up both sides of
	// the deletion point.
	// by the way, there is *always* a predecessor/successor to any collapsed
	// beach section, it's just impossible to have a collapsing first/last
	// beach sections on the beachline, since they obviously are unconstrained
	// on their left/right side.

	// look left
	lArc := previous.value.(*Beachsection)
	for lArc.circleEvent != nil &&
		abs_fn(x-lArc.circleEvent.x) < 1e-9 &&
		abs_fn(y-lArc.circleEvent.ycenter) < 1e-9 {

		previous = lArc.node.previous
		disappearingTransitions.appendLeft(lArc)
		s.detachBeachsection(lArc) // mark for reuse
		lArc = previous.value.(*Beachsection)
	}
	// even though it is not disappearing, I will also add the beach section
	// immediately to the left of the left-most collapsed beach section, for
	// convenience, since we need to refer to it later as this beach section
	// is the 'left' site of an edge for which a start point is set.
	disappearingTransitions.appendLeft(lArc)
	s.detachCircleEvent(lArc)

	// look right
	var rArc = next.value.(*Beachsection)
	for rArc.circleEvent != nil &&
		abs_fn(x-rArc.circleEvent.x) < 1e-9 &&
		abs_fn(y-rArc.circleEvent.ycenter) < 1e-9 {
		next = rArc.node.next
		disappearingTransitions.appendRight(rArc)
		s.detachBeachsection(rArc) // mark for reuse
		rArc = next.value.(*Beachsection)
	}
	// we also have to add the beach section immediately to the right of the
	// right-most collapsed beach section, since there is also a disappearing
	// transition representing an edge's start point on its left.
	disappearingTransitions.appendRight(rArc)
	s.detachCircleEvent(rArc)

	// walk through all the disappearing transitions between beach sections and
	// set the start point of their (implied) edge.
	nArcs := len(disappearingTransitions)

	for iArc := 1; iArc < nArcs; iArc++ {
		rArc = disappearingTransitions[iArc]
		lArc = disappearingTransitions[iArc-1]

		lSite := s.getCell(lArc.site)
		rSite := s.getCell(rArc.site)

		s.setEdgeStartpoint(rArc.edge, lSite, rSite, vertex)
	}

	// create a new edge as we have now a new transition between
	// two beach sections which were previously not adjacent.
	// since this edge appears as a new vertex is defined, the vertex
	// actually define an end point of the edge (relative to the site
	// on the left)
	lArc = disappearingTransitions[0]
	rArc = disappearingTransitions[nArcs-1]
	lSite := s.getCell(lArc.site)
	rSite := s.getCell(rArc.site)

	rArc.edge = s.createEdge(lSite, rSite, NoVertex, vertex)

	// create circle events if any for beach sections left in the beachline
	// adjacent to collapsed sections
	s.attachCircleEvent(lArc)
	s.attachCircleEvent(rArc)
}

func (s *Voronoi) addBeachsection(site Vertex) {
	x := site.X
	directrix := site.Y

	// find the left and right beach sections which will surround the newly
	// created beach section.
	// rhill 2011-06-01: This loop is one of the most often executed,
	// hence we expand in-place the comparison-against-epsilon calls.
	var lNode, rNode *rbNode
	var dxl, dxr float64
	node := s.beachline.root

	for node != nil {
		nodeBeachline := node.value.(*Beachsection)
		dxl = leftBreakPoint(nodeBeachline, directrix) - x
		// x lessThanWithEpsilon xl => falls somewhere before the left edge of the beachsection
		if dxl > 1e-9 {
			// this case should never happen
			// if (!node.rbLeft) {
			//    rNode = node.rbLeft;
			//    break;
			//    }
			node = node.left
		} else {
			dxr = x - rightBreakPoint(nodeBeachline, directrix)
			// x greaterThanWithEpsilon xr => falls somewhere after the right edge of the beachsection
			if dxr > 1e-9 {
				if node.right == nil {
					lNode = node
					break
				}
				node = node.right
			} else {
				// x equalWithEpsilon xl => falls exactly on the left edge of the beachsection
				if dxl > -1e-9 {
					lNode = node.previous
					rNode = node
				} else if dxr > -1e-9 {
					// x equalWithEpsilon xr => falls exactly on the right edge of the beachsection
					lNode = node
					rNode = node.next
					// falls exactly somewhere in the middle of the beachsection
				} else {
					lNode = node
					rNode = node
				}
				break
			}
		}
	}

	var lArc, rArc *Beachsection

	if lNode != nil {
		lArc = lNode.value.(*Beachsection)
	}
	if rNode != nil {
		rArc = rNode.value.(*Beachsection)
	}

	// at this point, keep in mind that lArc and/or rArc could be
	// undefined or null.

	// create a new beach section object for the site and add it to RB-tree
	newArc := &Beachsection{site: site}
	if lArc == nil {
		s.beachline.insertSuccessor(nil, newArc)
	} else {
		s.beachline.insertSuccessor(lArc.node, newArc)
	}

	// cases:
	//

	// [null,null]
	// least likely case: new beach section is the first beach section on the
	// beachline.
	// This case means:
	//   no new transition appears
	//   no collapsing beach section
	//   new beachsection become root of the RB-tree
	if lArc == nil && rArc == nil {
		return
	}

	// [lArc,rArc] where lArc == rArc
	// most likely case: new beach section split an existing beach
	// section.
	// This case means:
	//   one new transition appears
	//   the left and right beach section might be collapsing as a result
	//   two new nodes added to the RB-tree
	if lArc == rArc {
		// invalidate circle event of split beach section
		s.detachCircleEvent(lArc)

		// split the beach section into two separate beach sections
		rArc = &Beachsection{site: lArc.site}
		s.beachline.insertSuccessor(newArc.node, rArc)

		// since we have a new transition between two beach sections,
		// a new edge is born
		newArc.edge = s.createEdge(s.getCell(lArc.site), s.getCell(newArc.site), NoVertex, NoVertex)
		rArc.edge = newArc.edge

		// check whether the left and right beach sections are collapsing
		// and if so create circle events, to be notified when the point of
		// collapse is reached.
		s.attachCircleEvent(lArc)
		s.attachCircleEvent(rArc)
		return
	}

	// [lArc,null]
	// even less likely case: new beach section is the *last* beach section
	// on the beachline -- this can happen *only* if *all* the previous beach
	// sections currently on the beachline share the same y value as
	// the new beach section.
	// This case means:
	//   one new transition appears
	//   no collapsing beach section as a result
	//   new beach section become right-most node of the RB-tree
	if lArc != nil && rArc == nil {
		newArc.edge = s.createEdge(s.getCell(lArc.site), s.getCell(newArc.site), NoVertex, NoVertex)
		return
	}

	// [null,rArc]
	// impossible case: because sites are strictly processed from top to bottom,
	// and left to right, which guarantees that there will always be a beach section
	// on the left -- except of course when there are no beach section at all on
	// the beach line, which case was handled above.
	// rhill 2011-06-02: No point testing in non-debug version
	//if (!lArc && rArc) {
	//    throw "Voronoi.addBeachsection(): What is this I don't even";
	//    }

	// [lArc,rArc] where lArc != rArc
	// somewhat less likely case: new beach section falls *exactly* in between two
	// existing beach sections
	// This case means:
	//   one transition disappears
	//   two new transitions appear
	//   the left and right beach section might be collapsing as a result
	//   only one new node added to the RB-tree
	if lArc != rArc {
		// invalidate circle events of left and right sites
		s.detachCircleEvent(lArc)
		s.detachCircleEvent(rArc)

		// an existing transition disappears, meaning a vertex is defined at
		// the disappearance point.
		// since the disappearance is caused by the new beachsection, the
		// vertex is at the center of the circumscribed circle of the left,
		// new and right beachsections.
		// http://mathforum.org/library/drmath/view/55002.html
		// Except that I bring the origin at A to simplify
		// calculation
		LeftSite := lArc.site
		ax := LeftSite.X
		ay := LeftSite.Y
		bx := site.X - ax
		by := site.Y - ay
		RightSite := rArc.site
		cx := RightSite.X - ax
		cy := RightSite.Y - ay
		d := 2 * (bx*cy - by*cx)
		hb := bx*bx + by*by
		hc := cx*cx + cy*cy
		vertex := Vertex{(cy*hb-by*hc)/d + ax, (bx*hc-cx*hb)/d + ay}

		lCell := s.getCell(LeftSite)
		cell := s.getCell(site)
		rCell := s.getCell(RightSite)

		// one transition disappear
		s.setEdgeStartpoint(rArc.edge, lCell, rCell, vertex)

		// two new transitions appear at the new vertex location
		newArc.edge = s.createEdge(lCell, cell, NoVertex, vertex)
		rArc.edge = s.createEdge(cell, rCell, NoVertex, vertex)

		// check whether the left and right beach sections are collapsing
		// and if so create circle events, to handle the point of collapse.
		s.attachCircleEvent(lArc)
		s.attachCircleEvent(rArc)
		return
	}
}

type circleEvent struct {
	node    *rbNode
	site    Vertex
	arc     *Beachsection
	x       float64
	y       float64
	ycenter float64
}

func (s *circleEvent) bindToNode(node *rbNode) {
	s.node = node
}

func (s *circleEvent) getNode() *rbNode {
	return s.node
}

func (s *Voronoi) attachCircleEvent(arc *Beachsection) {
	lArc := arc.node.previous
	rArc := arc.node.next
	if lArc == nil || rArc == nil {
		return // does that ever happen?
	}
	LeftSite := lArc.value.(*Beachsection).site
	cSite := arc.site
	RightSite := rArc.value.(*Beachsection).site

	// If site of left beachsection is same as site of
	// right beachsection, there can't be convergence
	if LeftSite == RightSite {
		return
	}

	// Find the circumscribed circle for the three sites associated
	// with the beachsection triplet.
	// rhill 2011-05-26: It is more efficient to calculate in-place
	// rather than getting the resulting circumscribed circle from an
	// object returned by calling Voronoi.circumcircle()
	// http://mathforum.org/library/drmath/view/55002.html
	// Except that I bring the origin at cSite to simplify calculations.
	// The bottom-most part of the circumcircle is our Fortune 'circle
	// event', and its center is a vertex potentially part of the final
	// Voronoi diagram.
	bx := cSite.X
	by := cSite.Y
	ax := LeftSite.X - bx
	ay := LeftSite.Y - by
	cx := RightSite.X - bx
	cy := RightSite.Y - by

	// If points l->c->r are clockwise, then center beach section does not
	// collapse, hence it can't end up as a vertex (we reuse 'd' here, which
	// sign is reverse of the orientation, hence we reverse the test.
	// http://en.wikipedia.org/wiki/Curve_orientation#Orientation_of_a_simple_polygon
	// rhill 2011-05-21: Nasty finite precision error which caused circumcircle() to
	// return infinites: 1e-12 seems to fix the problem.
	d := 2 * (ax*cy - ay*cx)
	if d >= -2e-12 {
		return
	}

	ha := ax*ax + ay*ay
	hc := cx*cx + cy*cy
	x := (cy*ha - ay*hc) / d
	y := (ax*hc - cx*ha) / d
	ycenter := y + by

	// Important: ybottom should always be under or at sweep, so no need
	// to waste CPU cycles by checking

	// recycle circle event object if possible
	circleEventInst := &circleEvent{
		arc:     arc,
		site:    cSite,
		x:       x + bx,
		y:       ycenter + math.Sqrt(x*x+y*y),
		ycenter: ycenter,
	}

	arc.circleEvent = circleEventInst

	// find insertion point in RB-tree: circle events are ordered from
	// smallest to largest
	var predecessor *rbNode = nil
	node := s.circleEvents.root
	for node != nil {
		nodeValue := node.value.(*circleEvent)
		if circleEventInst.y < nodeValue.y || (circleEventInst.y == nodeValue.y && circleEventInst.x <= nodeValue.x) {
			if node.left != nil {
				node = node.left
			} else {
				predecessor = node.previous
				break
			}
		} else {
			if node.right != nil {
				node = node.right
			} else {
				predecessor = node
				break
			}
		}
	}
	s.circleEvents.insertSuccessor(predecessor, circleEventInst)
	if predecessor == nil {
		s.firstCircleEvent = circleEventInst
	}
}

func (s *Voronoi) detachCircleEvent(arc *Beachsection) {
	circle := arc.circleEvent
	if circle != nil {
		if circle.node.previous == nil {
			if circle.node.next != nil {
				s.firstCircleEvent = circle.node.next.value.(*circleEvent)
			} else {
				s.firstCircleEvent = nil
			}
		}
		s.circleEvents.removeNode(circle.node) // remove from RB-tree
		arc.circleEvent = nil
	}
}

// Bounding Box
type BBox struct {
	Xl,
	Yt,
	Xr,
	Yb float64
}

// Create new Bounding Box by providing top-left and bottom-right coordinates
func NewBBox(xLeft, yTop, xRight, yBottom float64) BBox {
	return BBox{xLeft, yTop, xRight, yBottom}
}

// connect dangling edges (not if a cursory test tells us
// it is not going to be visible.
// return value:
//   false: the dangling endpoint couldn't be connected
//   true: the dangling endpoint could be connected
func connectEdge(edge *Edge, bbox BBox) bool {
	// skip if end point already connected
	vb := edge.Vb.Vertex
	if vb != NoVertex {
		return true
	}

	// make local copy for performance purpose
	va := edge.Va.Vertex
	xl := bbox.Xl
	xr := bbox.Xr
	yt := bbox.Yt
	yb := bbox.Yb
	LeftSite := edge.LeftCell.Site
	RightSite := edge.RightCell.Site
	lx := LeftSite.X
	ly := LeftSite.Y
	rx := RightSite.X
	ry := RightSite.Y
	fx := (lx + rx) / 2
	fy := (ly + ry) / 2

	var fm, fb float64

	// if we reach here, this means cells which use this edge will need
	// to be closed, whether because the edge was removed, or because it
	// was connected to the bounding box.
	// TODO(zzwx): Is it what it is?
	edge.LeftCell.closeMe = true
	edge.RightCell.closeMe = true
	//this.cells[lSite.voronoiId].closeMe = true;
	//this.cells[rSite.voronoiId].closeMe = true;

	// get the line equation of the bisector if line is not vertical
	if !equalWithEpsilon(ry, ly) {
		fm = (lx - rx) / (ry - ly)
		fb = fy - fm*fx
	}

	// remember, direction of line (relative to left site):
	// upward: left.X < right.X
	// downward: left.X > right.X
	// horizontal: left.X == right.X
	// upward: left.X < right.X
	// rightward: left.Y < right.Y
	// leftward: left.Y > right.Y
	// vertical: left.Y == right.Y

	// depending on the direction, find the best side of the
	// bounding box to use to determine a reasonable start point

	// special case: vertical line
	if equalWithEpsilon(ry, ly) {
		// doesn't intersect with viewport
		if fx < xl || fx >= xr {
			return false
		}
		// downward
		if lx > rx {
			if va == NoVertex {
				va = Vertex{fx, yt}
			} else if va.Y >= yb {
				return false
			}
			vb = Vertex{fx, yb}
			// upward
		} else {
			if va == NoVertex {
				va = Vertex{fx, yb}
			} else if va.Y < yt {
				return false
			}
			vb = Vertex{fx, yt}
		}
		// closer to vertical than horizontal, connect start point to the
		// top or bottom side of the bounding box
	} else if fm < -1 || fm > 1 {
		// downward
		if lx > rx {
			if va == NoVertex {
				va = Vertex{(yt - fb) / fm, yt}
			} else if va.Y >= yb {
				return false
			}
			vb = Vertex{(yb - fb) / fm, yb}
			// upward
		} else {
			if va == NoVertex {
				va = Vertex{(yb - fb) / fm, yb}
			} else if va.Y < yt {
				return false
			}
			vb = Vertex{(yt - fb) / fm, yt}
		}
		// closer to horizontal than vertical, connect start point to the
		// left or right side of the bounding box
	} else {
		// rightward
		if ly < ry {
			if va == NoVertex {
				va = Vertex{xl, fm*xl + fb}
			} else if va.X >= xr {
				return false
			}
			vb = Vertex{xr, fm*xr + fb}
			// leftward
		} else {
			if va == NoVertex {
				va = Vertex{xr, fm*xr + fb}
			} else if va.X < xl {
				return false
			}
			vb = Vertex{xl, fm*xl + fb}
		}
	}
	edge.Va.Vertex = va
	edge.Vb.Vertex = vb
	return true
}

// line-clipping code taken from:
//   Liang-Barsky function by Daniel White
//   http://www.skytopia.com/project/articles/compsci/clipping.html
// Thanks!
// A bit modified to minimize code paths
func clipEdge(edge *Edge, bbox BBox) bool {
	ax := edge.Va.X
	ay := edge.Va.Y
	bx := edge.Vb.X
	by := edge.Vb.Y
	t0 := float64(0)
	t1 := float64(1)
	dx := bx - ax
	dy := by - ay

	// left
	q := ax - bbox.Xl
	if dx == 0 && q < 0 {
		return false
	}
	r := -q / dx
	if dx < 0 {
		if r < t0 {
			return false
		}
		if r < t1 {
			t1 = r
		}
	} else if dx > 0 {
		if r > t1 {
			return false
		}
		if r > t0 {
			t0 = r
		}
	}
	// right
	q = bbox.Xr - ax
	if dx == 0 && q < 0 {
		return false
	}
	r = q / dx
	if dx < 0 {
		if r > t1 {
			return false
		}
		if r > t0 {
			t0 = r
		}
	} else if dx > 0 {
		if r < t0 {
			return false
		}
		if r < t1 {
			t1 = r
		}
	}

	// top
	q = ay - bbox.Yt
	if dy == 0 && q < 0 {
		return false
	}
	r = -q / dy
	if dy < 0 {
		if r < t0 {
			return false
		}
		if r < t1 {
			t1 = r
		}
	} else if dy > 0 {
		if r > t1 {
			return false
		}
		if r > t0 {
			t0 = r
		}
	}
	// bottom
	q = bbox.Yb - ay
	if dy == 0 && q < 0 {
		return false
	}
	r = q / dy
	if dy < 0 {
		if r > t1 {
			return false
		}
		if r > t0 {
			t0 = r
		}
	} else if dy > 0 {
		if r < t0 {
			return false
		}
		if r < t1 {
			t1 = r
		}
	}

	// if we reach this point, Voronoi edge is within bbox

	// if t0 > 0, va needs to change
	// rhill 2011-06-03: we need to create a new vertex rather
	// than modifying the existing one, since the existing
	// one is likely shared with at least another edge
	if t0 > 0 {
		edge.Va.Vertex = Vertex{ax + t0*dx, ay + t0*dy}
	}

	// if t1 < 1, vb needs to change
	// rhill 2011-06-03: we need to create a new vertex rather
	// than modifying the existing one, since the existing
	// one is likely shared with at least another edge
	if t1 < 1 {
		edge.Vb.Vertex = Vertex{ax + t1*dx, ay + t1*dy}
	}

	// va and/or vb were clipped, thus we will need to close
	// cells which use this edge.
	if t0 > 0 || t1 < 1 {
		// TODO(zzwx): Is this what it is?
		edge.LeftCell.closeMe = true
		edge.RightCell.closeMe = true
		//this.cells[edge.lSite.voronoiId].closeMe = true;
		//this.cells[edge.rSite.voronoiId].closeMe = true;
	}

	return true
}

func equalWithEpsilon(a, b float64) bool {
	return math.Abs(a-b) < 1e-9
}

func lessThanWithEpsilon(a, b float64) bool {
	return b-a > 1e-9
}

func greaterThanWithEpsilon(a, b float64) bool {
	return a-b > 1e-9
}

// Connect/cut edges at bounding box
func (s *Voronoi) clipEdges(bbox BBox) {
	// connect all dangling edges to bounding box
	// or get rid of them if it can't be done
	abs_fn := math.Abs

	// iterate backward so we can splice safely
	for iEdge := len(s.edges) - 1; iEdge >= 0; iEdge-- {
		edge := s.edges[iEdge]
		// edge is removed if:
		//   it is wholly outside the bounding box
		//   it is looking more like a point than a line
		if !connectEdge(edge, bbox) ||
			!clipEdge(edge, bbox) ||
			(abs_fn(edge.Va.X-edge.Vb.X) < 1e-9 && abs_fn(edge.Va.Y-edge.Vb.Y) < 1e-9) {
			edge.Va.Vertex = NoVertex
			edge.Vb.Vertex = NoVertex
			SpliceEdges(&s.edges, iEdge, 1)
			//s.edges[iEdge] = s.edges[len(s.edges)-1]
			//s.edges = s.edges[0 : len(s.edges)-1]
		}
	}
}

// closeCells closes the cells.
// The cells are bound by the supplied bounding box.
// Each cell refers to its associated site, and a list
// of halfedges ordered counterclockwise.
func (s *Voronoi) closeCells(bbox BBox) {
	xl := bbox.Xl
	xr := bbox.Xr
	yt := bbox.Yt
	yb := bbox.Yb
	cells := s.cells
	abs_fn := math.Abs

	for iCell := len(cells) - 1; iCell >= 0; iCell-- {
		cell := cells[iCell]
		// prune, order halfedges counterclockwise, then add missing ones
		// required to close cells
		if cell.prepareHalfedges() == 0 {
			continue
		}
		if !cell.closeMe {
			continue
		}
		// find first 'unclosed' point.
		// an 'unclosed' point will be the end point of a halfedge which
		// does not match the start point of the following halfedge
		//halfedges := cell.Halfedges
		//nHalfedges := len(cell.Halfedges)
		// special case: only one site, in which case, the viewport is the cell
		// ...

		// all other cases
		iLeft := 0
		for iLeft < len(cell.Halfedges) {
			va := cell.Halfedges[iLeft].GetEndpoint()
			vz := cell.Halfedges[(iLeft+1)%len(cell.Halfedges)].GetStartpoint()
			// if end point is not equal to start point, we need to add the missing
			// halfedge(s) to close the cell
			if abs_fn(va.X-vz.X) >= 1e-9 || abs_fn(va.Y-vz.Y) >= 1e-9 {

				// rhill 2013-12-02:
				// "Holes" in the halfedges are not necessarily always adjacent.
				// https://github.com/gorhill/Javascript-Voronoi/issues/16

				// find entry point:
				switch {
				// walk downward along left side
				case equalWithEpsilon(va.X, xl) && lessThanWithEpsilon(va.Y, yb):
					lastBorderSegment := equalWithEpsilon(vz.X, xl)
					var vb Vertex
					if lastBorderSegment {
						vb = Vertex{xl, vz.Y}
					} else {
						vb = Vertex{xl, yb}
					}
					edge := s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					va = vb
					fallthrough
				// walk rightward along bottom side
				case equalWithEpsilon(va.Y, yb) && lessThanWithEpsilon(va.X, xr):
					lastBorderSegment := equalWithEpsilon(vz.Y, yb)
					var vb Vertex
					if lastBorderSegment {
						vb = Vertex{vz.X, yb}
					} else {
						vb = Vertex{xr, yb}
					}
					edge := s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					va = vb
					fallthrough
				// walk upward along right side
				case equalWithEpsilon(va.X, xr) && greaterThanWithEpsilon(va.Y, yt):
					lastBorderSegment := equalWithEpsilon(vz.X, xr)
					var vb Vertex
					if lastBorderSegment {
						vb = Vertex{xr, vz.Y}
					} else {
						vb = Vertex{xr, yt}
					}
					edge := s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					va = vb
					fallthrough
				// walk leftward along top side
				case equalWithEpsilon(va.Y, yt) && greaterThanWithEpsilon(va.X, xl):
					lastBorderSegment := equalWithEpsilon(vz.Y, yt)
					var vb Vertex
					if lastBorderSegment {
						vb = Vertex{vz.X, yt}
					} else {
						vb = Vertex{xl, yt}
					}
					edge := s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					va = vb

					// walk downward along left side
					lastBorderSegment = equalWithEpsilon(vz.X, xl)
					if lastBorderSegment {
						vb = Vertex{xl, vz.Y}
					} else {
						vb = Vertex{xl, yb}
					}
					edge = s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					va = vb

					// walk rightward along bottom side
					lastBorderSegment = equalWithEpsilon(vz.Y, yb)
					if lastBorderSegment {
						vb = Vertex{vz.X, yb}
					} else {
						vb = Vertex{xr, yb}
					}
					edge = s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					va = vb

					// walk upward along right side
					lastBorderSegment = equalWithEpsilon(vz.X, xr)
					if lastBorderSegment {
						vb = Vertex{xr, vz.Y}
					} else {
						vb = Vertex{xr, yt}
					}
					edge = s.createBorderEdge(cell, va, vb)
					iLeft++
					SpliceHalfedges(&cell.Halfedges, iLeft, 0, newHalfedge(edge, cell, nil))
					//nHalfedges++
					if lastBorderSegment {
						break
					}
					// no va = vb

					fallthrough
				default:
					// TODO(zzwx): graceful
					panic("Voronoi.closeCells() > this makes no sense!")
				}
			}
			iLeft++
		}
		cell.closeMe = false
	}
}

func (s *Voronoi) gatherVertexEdges() {
	vertexEdgeMap := make(map[Vertex][]*Edge)

	for _, edge := range s.edges {
		vertexEdgeMap[edge.Va.Vertex] = append(
			vertexEdgeMap[edge.Va.Vertex], edge)
		vertexEdgeMap[edge.Vb.Vertex] = append(
			vertexEdgeMap[edge.Vb.Vertex], edge)
	}

	for vertex, edgeSlice := range vertexEdgeMap {
		for _, edge := range edgeSlice {
			if vertex == edge.Va.Vertex {
				edge.Va.Edges = edgeSlice
			}
			if vertex == edge.Vb.Vertex {
				edge.Vb.Edges = edgeSlice
			}
		}
	}
}

// Compute voronoi diagram. If closeCells == true, edges from bounding box will be
// included in diagram.
func ComputeDiagram(sites []Vertex, bbox BBox, closeCells bool) *Diagram {
	s := &Voronoi{
		cellsMap: make(map[Vertex]*Cell),
	}

	siteEvents := make([]Vertex, len(sites))
	for i := 0; i < len(sites); i++ {
		siteEvents[i] = sites[i]
	}

	// Initialize site event queue
	sort.Sort(VerticesByYX{siteEvents})

	pop := func() *Vertex {
		if len(siteEvents) == 0 {
			return nil
		}
		site := siteEvents[len(siteEvents)-1]
		siteEvents = siteEvents[:len(siteEvents)-1]
		return &site
	}

	// process queue
	site := pop()
	invalid := true // emulate undefined: first xsitex, xsitey are invalid. To address https://github.com/gorhill/Javascript-Voronoi/pull/14
	var xsitex float64
	var xsitey float64
	var circle *circleEvent

	// main loop
	for {
		// we need to figure whether we handle a site or circle event
		// for this we find out if there is a site event and it is
		// 'earlier' than the circle event
		circle = s.firstCircleEvent

		// add beach section
		if site != nil && (circle == nil || site.Y < circle.y || (site.Y == circle.y && site.X < circle.x)) {
			// only if site is not a duplicate
			if invalid || site.X != xsitex || site.Y != xsitey {
				// first create cell for new site
				nCell := newCell(*site)
				s.cells = append(s.cells, nCell)
				s.cellsMap[*site] = nCell
				// then create a beachsection for that site
				s.addBeachsection(*site)
				// remember last site coords to detect duplicate
				xsitey = site.Y
				xsitex = site.X
				invalid = false
			}
			site = pop()
			// remove beach section
		} else if circle != nil {
			s.removeBeachsection(circle.arc)
			// all done, quit
		} else {
			break
		}
	}

	// wrapping-up:
	//   connect dangling edges to bounding box
	//   cut edges as per bounding box
	//   discard edges completely outside bounding box
	//   discard edges which are point-like
	s.clipEdges(bbox)

	//   add missing edges in order to close opened cells
	if closeCells {
		s.closeCells(bbox)
	} else {
		for _, cell := range s.cells {
			cell.prepareHalfedges()
		}
	}

	// TODO(zzwx): Check if this is still necessary
	s.gatherVertexEdges()

	result := &Diagram{
		Edges: s.edges,
		Cells: s.cells,
	}
	return result
}
