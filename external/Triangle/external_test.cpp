#ifdef SINGLE
#define REAL float
#else /* not SINGLE */
#define REAL double
#endif /* not SINGLE */

#include <math.h>
#include <iostream>
#include <vector>

/* The vertex data structure.  Each vertex is actually an array of REALs.    */
/*   The number of REALs is unknown until runtime.  An integer boundary      */
/*   marker, and sometimes a pointer to a triangle, is appended after the    */
/*   REALs.                                                                  */
typedef REAL *vertex;


struct Polygon {
    int nVertices;
    REAL* vertices;
    REAL res;

    Polygon(int nVertices, REAL* vertices, REAL res) :
    nVertices(nVertices),
    vertices(vertices),
    res(res)
    {}
};


bool pointInPolygon(REAL x, REAL y, int nVertices, REAL* vertices) {
    // Polygon must be defined in an oriented fashion. See https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html
    int j = nVertices - 1;
    bool inside = false;
    bool flippedYet = false;

    for (int i = 0; i < nVertices; i++) {
        if ((vertices[2 * i + 1] > y) != (vertices[2 * j + 1] > y)) {
            REAL xx = (vertices[2 * j] - vertices[2 * i]) * (y - vertices[2 * i + 1]) / (vertices[2 * j + 1] - vertices[2 * i + 1]) + vertices[2 * i];
            if (x < xx) {
                inside = !inside;
            }
        }
        j = i;
    }

    return inside;
}


// By default, return a maximum "resolution" that will always be large.
REAL defaultFunction(REAL x, REAL y, REAL* args) {
    return INFINITY;
}


struct Callable {
    std::vector<REAL*> arg_list;
    std::vector<REAL (*)(REAL, REAL, REAL*)> function_list;
    std::vector<Polygon> polygon_list;

    REAL* args; // default arguments, never used

    REAL operator()(REAL x, REAL y) {
        // Return the default function if there are no constraints defined
        if (function_list.size() == 0 && polygon_list.size() == 0) {
            return defaultFunction(x, y, args);
        }

        REAL maxres = INFINITY;

        // Loop through the defined functions and take their minimum envelope
        for (int i = 0; i < function_list.size(); i++) {
            REAL testres = function_list[i](x, y, arg_list[i]);
            if (testres < maxres) { maxres = testres; }
        }

        // Loop through the defined polygons and set (minres, maxres) for
        // (interior, exterior)
        for (int i = 0; i < polygon_list.size(); i++) {
            if (pointInPolygon(x, y, polygon_list[i].nVertices, polygon_list[i].vertices)) {
                if (polygon_list[i].res < maxres) {
                    maxres = polygon_list[i].res;
                }
            }
        }

        return maxres;
    }
} callable;


extern "C" int triunsuitable(vertex triorg, vertex tridest, vertex triapex, REAL area)
{
    REAL dxoa, dxda, dxod;
    REAL dyoa, dyda, dyod;
    REAL oalen, dalen, odlen;
    REAL maxlen;

    dxoa = triorg[0] - triapex[0];
    dyoa = triorg[1] - triapex[1];
    dxda = tridest[0] - triapex[0];
    dyda = tridest[1] - triapex[1];
    dxod = triorg[0] - tridest[0];
    dyod = triorg[1] - tridest[1];
    /* Find the squares of the lengths of the triangle's three edges. */
    oalen = dxoa * dxoa + dyoa * dyoa;
    dalen = dxda * dxda + dyda * dyda;
    odlen = dxod * dxod + dyod * dyod;
    /* Find the square of the length of the longest edge. */
    maxlen = (dalen > oalen) ? dalen : oalen;
    maxlen = (odlen > maxlen) ? odlen : maxlen;

    REAL trix = (triorg[0] + tridest[0] + triapex[0]) / 3.0;
    REAL triy = (triorg[1] + tridest[1] + triapex[1]) / 3.0;
    REAL maxres = callable(trix, triy);

    if (maxlen > maxres * maxres) { return 1; }
    else { return 0; }
}


extern "C" void setTriangleResolutionFunction(REAL (*f)(REAL, REAL, REAL*), REAL* args) {
    callable.function_list.push_back(f);
    callable.arg_list.push_back(args);
    return;
}


extern "C" void setTriangleResolutionIndicatorFunction(int nVertices, REAL* vertices, REAL res) {
    Polygon polygon(nVertices, vertices, res);
    callable.polygon_list.push_back(polygon);
    return;
}


extern "C" void clearTriangleResolutionFunction() {
    callable.function_list.clear();
    callable.arg_list.clear();
    return;
}
