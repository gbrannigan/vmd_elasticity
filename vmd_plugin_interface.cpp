#include <tcl.h>
#include "bending_modulus.h"

extern "C" {
    int Vmd_bending_modulus_Init(Tcl_Interp *interp);
}
static int obj_calc_bending_modulus(ClientData, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[]);

int Vmd_bending_modulus_Init(Tcl_Interp *interp) {
    Tcl_CreateObjCommand(interp, "calc_bending_modulus", obj_calc_bending_modulus, (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL);
    Tcl_EvalEx(interp, "package provide bending_modulus 1.0", -1, 0);
    return TCL_OK;
}

int tcl_printf(Tcl_Interp *interp, const char *fmt, ...) {
    char buf[8192];
    va_list va;
    va_start(va, fmt);
    vsnprintf(buf, 8192, fmt, va);
    va_end(va);
    Tcl_Obj *msg = Tcl_ObjPrintf ("puts \"bending_modulus: %s\"", buf);
    int result = Tcl_EvalObjEx(interp, msg, TCL_EVAL_DIRECT);
    if (result != TCL_OK) { return TCL_ERROR; }
}

// Parse a vector of floats from a Tcl object
// This function flagrantly stolen from Jerome Henin's qwrap
int parse_vector (Tcl_Obj * const obj, std::vector<float> &vec, Tcl_Interp *interp)
{
    Tcl_Obj **data;
    int num;
    double d;

    if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "bending_modulus: parse_vector: error parsing arguments", TCL_STATIC);
        return -1;
    }

    vec.resize(num);

    for (int i = 0; i < num; i++) {
        if (Tcl_GetDoubleFromObj(interp, data[i], &d) != TCL_OK) {
            Tcl_SetResult(interp, (char *) "bending_modulus: parse_vector: error parsing vector element as floating-point", TCL_STATIC);
            return -1;
        }
        // Tcl gives us doubles, make them float
        vec[i] = float (d);
    }
    return num;
}

// Parse a vector of integers from a Tcl object
// Also flagrantly stolen from Jerome Henin
int parse_ivector (Tcl_Obj * const obj, std::vector<int> &vec, Tcl_Interp *interp, bool fromDouble)
{
    Tcl_Obj **data;
    int num;
    double d;

    if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "bending_modulus: parse_ivector: error parsing arguments", TCL_STATIC);
        return -1;
    }

    vec.resize(num);

    if (fromDouble == false) {
        for (int i = 0; i < num; i++) {
            if (Tcl_GetIntFromObj(interp, data[i], &vec[i]) != TCL_OK) {
                Tcl_SetResult(interp, (char *) "bending_modulus: parse_ivector: error parsing vector element as integer", TCL_STATIC);
                return -1;
            }
        }
    } else {
        // do a double-to-int conversion first
        for (int i = 0; i < num; i++) {
            if (Tcl_GetDoubleFromObj(interp, data[i], &d) != TCL_OK) {
                Tcl_SetResult(interp, (char *) "bending_modulus: parse_ivector: error parsing vector element as integer", TCL_STATIC);
                return -1;
            }
            vec[i] = int (d);
        }
    }
    return num;
}

// Retrieves atom indices of an existing atomselect object
int get_atomselect_indices(Tcl_Interp *interp, Tcl_Obj *atomselect, std::vector<int> &out)
{
    Tcl_Obj *script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendToObj(script, " get index", -1);
    int result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
    if(result != TCL_OK)
    {
        Tcl_SetResult(interp, (char *) "bending_modulus: Error getting atom indices from atomselect object", TCL_STATIC);
        return -1;
    }
    return parse_ivector(Tcl_GetObjResult(interp), out, interp, false);
}


// The actual bending_modulus command
static int obj_calc_bending_modulus(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[]) {
    tcl_printf(interp, "Hello, VMD!");
    Tcl_SetResult(interp, (char *) "calc_bending_modulus: There was some kind of error. Please check into this.",
                  TCL_STATIC);
    return TCL_ERROR;

}