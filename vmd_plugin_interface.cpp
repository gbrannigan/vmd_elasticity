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
int parse_vector(Tcl_Obj * const obj, std::vector<float> &vec, Tcl_Interp *interp) {

    Tcl_Obj **data;
    int num;
    if (Tcl_ListObjGetElements(interp, obj, &num, &data) != TCL_OK) {
        Tcl_SetResult(interp, (char *) "bending_modulus: parse_vector: error parsing arguments", TCL_STATIC);
        return -1;
    }

    vec.resize(num);

    for (int i = 0; i < num; i++) {
        double d;
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
int parse_ivector(Tcl_Obj * const obj, std::vector<int> &vec, Tcl_Interp *interp, bool fromDouble) {
    Tcl_Obj **data;
    int num;

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
            double d;
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
int get_atomselect_indices(Tcl_Interp *interp, Tcl_Obj *atomselect, std::vector<int> &out) {
    Tcl_Obj *script = Tcl_DuplicateObj(atomselect);
    Tcl_AppendToObj(script, " get index", -1);
    int result = Tcl_EvalObjEx(interp, script, TCL_EVAL_DIRECT);
    if(result != TCL_OK) {
        Tcl_SetResult(interp, (char *) "bending_modulus: Error getting atom indices from atomselect object", TCL_STATIC);
        return -1;
    }
    return parse_ivector(Tcl_GetObjResult(interp), out, interp, false);
}

// Runs a TCL command that returns a Tcl object. Returns that object, and the status code in &status
Tcl_Obj *run_tcl_command_returning_obj(Tcl_Interp *interp, const char *text, int &status) {
    Tcl_Obj *cmd =  Tcl_ObjPrintf(text);
    status = Tcl_EvalObjEx(interp, cmd, TCL_EVAL_DIRECT);
    if (status != TCL_OK) {
        Tcl_SetResult(interp, (char *) "bending_modulus: Something went wrong. This is an awful error message.", TCL_STATIC);
        return nullptr;
    }
    return Tcl_GetObjResult(interp);
}

// Runs a TCL command that returns a Tcl object. Returns that object, and the status code in &status
int run_tcl_command_returning_int(Tcl_Interp *interp, const char *text, int &status) {
    Tcl_Obj *obj = run_tcl_command_returning_obj(interp, text, status);
    if(status != TCL_OK) {
        return -1;
    }
    int data;
    status = Tcl_GetIntFromObj(interp, obj, &data);
    if (status != TCL_OK) {
        tcl_printf(interp, "hello %d. (TCL_OK = %d, TCL_ERROR = %d)", status, TCL_OK, TCL_ERROR);
        // Tcl_SetResult(interp, (char *) "bending_modulus: error parsing number of frames", TCL_STATIC);
        return -1;
    }
    return data;
}


// The actual bending_modulus command
static int obj_calc_bending_modulus(ClientData data, Tcl_Interp *interp, int argc, Tcl_Obj * const objv[]) {
    // TODO: Take an atomselect object containing only heads and tails
    // FIXME: For now, just uses the "atomselect top"

    // Get number of frames in the system
    int status;
    int num_frames = run_tcl_command_returning_int(interp, "molinfo top get numframes", status);
    tcl_printf(interp, "There are %d frames.", num_frames);

    FloatArray2D box(num_frames, 3);

    Tcl_Obj *atomselect_heads = run_tcl_command_returning_obj(interp, "atomselect top \"lipid and name P\"", status);
    Tcl_IncrRefCount(atomselect_heads);
    Tcl_Obj *atomselect_tails = run_tcl_command_returning_obj(interp, "atomselect top \"lipid and name C218\"", status);
    Tcl_IncrRefCount(atomselect_tails);

    std::vector<int> heads_indices, tails_indices;
    int num_heads = get_atomselect_indices(interp, atomselect_heads, heads_indices);
    int num_tails = get_atomselect_indices(interp, atomselect_tails, tails_indices);

    if(num_heads != num_tails) {
        Tcl_SetResult(interp, (char *) "bending_modulus: Number of heads does not equal number of tails", TCL_STATIC);
        return TCL_ERROR;
    }

    // Lipid coordinates array is arranged:
    //  0: head_x, head_y, head_z
    //  1: tail_x, tail_y, tail_z
    //  2: head_x, head_y, head_z
    //  3: tail_x, tail_y, tail_z
    //  ...
    FloatArray2D lipid(num_frames*num_heads*2, 3);

    for(int frame = 0; frame < num_frames; frame++) {
        // I guess we have to change the frame like this?
        Tcl_Obj *chgframe = Tcl_DuplicateObj(atomselect_heads);
        Tcl_AppendPrintfToObj (chgframe, " frame %i", frame);
        status = Tcl_EvalObjEx(interp, chgframe, TCL_EVAL_DIRECT);
        if(status != TCL_OK) {
            return TCL_ERROR;
        }

        chgframe = Tcl_DuplicateObj(atomselect_tails);
        Tcl_AppendPrintfToObj (chgframe, " frame %i", frame);
        status = Tcl_EvalObjEx(interp, chgframe, TCL_EVAL_DIRECT);
        if(status != TCL_OK) {
            return TCL_ERROR;
        }

        // Choose the molinfo in the frame we want. I guess we do this to get box info
        Tcl_Obj *mol_chgframe = Tcl_ObjPrintf ("molinfo top set frame %i", frame);
        status = Tcl_EvalObjEx(interp, mol_chgframe, TCL_EVAL_DIRECT);
        if(status != TCL_OK) {
            return TCL_ERROR;
        }

        // Actually get the box size
        Tcl_Obj *get_abc = Tcl_ObjPrintf ("molinfo top get {a b c}");
        status = Tcl_EvalObjEx(interp, get_abc, TCL_EVAL_DIRECT);
        if(status != TCL_OK) {
            return TCL_ERROR;
        }
        Tcl_Obj *object = Tcl_GetObjResult(interp);
        // Extract those three precious floats from the Tcl object
        std::vector<float> PBC;
        int num = parse_vector(object, PBC, interp);
        // Ensure there are three numbers, and they are orthogonal
        if (num != 3 || PBC[0]*PBC[1]*PBC[2] == 0.0) {
            Tcl_SetResult(interp, (char *) "bending_modulus: PBC is either not 3 numbers or not orthogonal", TCL_STATIC);
            return TCL_ERROR;
        }
        box.set(frame, 0, PBC[0]);
        box.set(frame, 1, PBC[1]);
        box.set(frame, 2, PBC[2]);

        // Get binary array of coordinates for this frame. I guess we get all of them, but we have the
        // indices of the heads and tails from above, so we can just save the ones we want.
        Tcl_Obj *get_ts = Tcl_ObjPrintf ("gettimestep %s %i", "top", frame);
        status = Tcl_EvalObjEx(interp, get_ts, TCL_EVAL_DIRECT);
        if (status != TCL_OK) {
            Tcl_SetResult(interp, (char *) "bending_modulus: error getting coordinates", TCL_STATIC);
            return TCL_ERROR;
        }

        Tcl_Obj *bytes = Tcl_GetObjResult(interp);
        Tcl_IncrRefCount(bytes);
        Tcl_InvalidateStringRep (bytes);
        int length;
        // cerr << "DEBUG bending_modulus: Got " << length << " bytes" << endl;
        float *coords = reinterpret_cast<float *> (Tcl_GetByteArrayFromObj(bytes, &length));

        // Seems like the coords are arranged like: xyzxyzxyz...
        // For each head/tail index, copy the x, y, z coords individually
        // Each lipid[xyz] array is stored as hthththt...
        int this_frame_offset = frame * num_heads * 2;
        for(int i = 0; i < num_heads; i++) {
            int idx = heads_indices[i];

            lipid.set(this_frame_offset + 2*i, 0, coords[idx*3]);
            lipid.set(this_frame_offset + 2*i, 1, coords[idx*3+1]);
            lipid.set(this_frame_offset + 2*i, 2, coords[idx*3+2]);

            idx = tails_indices[i];
            lipid.set(this_frame_offset + 2*i+1, 0, coords[idx*3]);
            lipid.set(this_frame_offset + 2*i+1, 1, coords[idx*3+1]);
            lipid.set(this_frame_offset + 2*i+1, 2, coords[idx*3+2]);
        }

    }

    // TODO: Do we need to decrement the tcl reference count?
    // Tcl_DecrRefCount(atomselect_heads);
    // Tcl_DecrRefCount(atomselect_tails);
    tcl_printf(interp, "Done loading frames. See stdout for output (for now).");

    // Fasten your seatbelts
    do_bending_modulus(lipid, box, num_heads, 12);

    Tcl_SetResult(interp, (char *) "bending_modulus: Implementation not finished yet.",
                  TCL_STATIC);
    return TCL_ERROR;

}