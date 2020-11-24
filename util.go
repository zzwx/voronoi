package voronoi

// Copied from github.com/zzwx/splice and updated for necessary types

func SpliceHalfedges(source *[]*Halfedge, start int, deleteCount int, item ...*Halfedge) (arrDeletedItems []*Halfedge) {
	arrDeletedItems = []*Halfedge{}
	if start > len(*source) {
		start = len(*source)
	}
	if start < 0 {
		start = len(*source) + start
	}
	if start < 0 {
		start = 0
	}
	if deleteCount < 0 {
		deleteCount = 0
	}
	if deleteCount > 0 {
		for i := 0; i < deleteCount; i++ {
			if i+start < len(*source) {
				arrDeletedItems = append(arrDeletedItems, (*source)[i+start])
			}
		}
	}
	deleteCount = len(arrDeletedItems) // Adjust to actual delete count
	grow := len(item) - deleteCount
	switch {
	case grow > 0: // So we grow
		*source = append(*source, make([]*Halfedge, grow)...)
		copy((*source)[start+deleteCount+grow:], (*source)[start+deleteCount:])
	case grow < 0: // So we shrink
		from := start + len(item)
		to := start + deleteCount
		copy((*source)[from:], (*source)[to:])
		*source = (*source)[:len(*source)+grow]
	}
	copy((*source)[start:], item)
	return
}

func SpliceEdges(source *[]*Edge, start int, deleteCount int, item ...*Edge) (arrDeletedItems []*Edge) {
	arrDeletedItems = []*Edge{}
	if start > len(*source) {
		start = len(*source)
	}
	if start < 0 {
		start = len(*source) + start
	}
	if start < 0 {
		start = 0
	}
	if deleteCount < 0 {
		deleteCount = 0
	}
	if deleteCount > 0 {
		for i := 0; i < deleteCount; i++ {
			if i+start < len(*source) {
				arrDeletedItems = append(arrDeletedItems, (*source)[i+start])
			}
		}
	}
	deleteCount = len(arrDeletedItems) // Adjust to actual delete count
	grow := len(item) - deleteCount
	switch {
	case grow > 0: // So we grow
		*source = append(*source, make([]*Edge, grow)...)
		copy((*source)[start+deleteCount+grow:], (*source)[start+deleteCount:])
	case grow < 0: // So we shrink
		from := start + len(item)
		to := start + deleteCount
		copy((*source)[from:], (*source)[to:])
		*source = (*source)[:len(*source)+grow]
	}
	copy((*source)[start:], item)
	return
}
