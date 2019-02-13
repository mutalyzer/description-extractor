#include <Python.h>

#include "extractor-core/include/extractor.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>


// Helper function to build a point position object
static PyObject const*
point_location(size_t const point)
{
    return Py_BuildValue("{s:s,s:i}", "type", "point", "position", point);
} // point_location


// Build a range location object. All variants are described using ranges.
static PyObject const*
range_location(size_t const start, size_t const end)
{
    PyObject const* const start_object = point_location(start);
    if (start_object == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if
    PyObject const* const end_object = point_location(end);
    if (end_object == NULL)
    {
        Py_DECREF(start_object);
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if
    return Py_BuildValue("{s:s,s:O,s:O}", "type", "range", "start", start_object, "end", end_object);
} // range_location


// Build a inserted object. This object holds the inserted part of a variant.
static PyObject const*
insertion_dict(mutalyzer::Variant const &variant)
{

    // Transposition (taken from the reference somewhere)
    if ((variant.type & mutalyzer::IDENTITY) == mutalyzer::IDENTITY)
    {
        PyObject const* const range = range_location(variant.transposition_start, variant.transposition_end);
        if (range == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        return Py_BuildValue("{s:s,s:O}", "source", "reference", "location", range);
    } // if

    // Inverse transposition (the reverse complement taken from the reference somewhere)
    else if ((variant.type & mutalyzer::REVERSE_COMPLEMENT) == mutalyzer::REVERSE_COMPLEMENT)
    {
        PyObject const* const range = range_location(variant.transposition_start, variant.transposition_end);
        if (range == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        return Py_BuildValue("{s:b,s:s,s:O}", "inverted", true, "source", "reference", "location", range);
    } // if

    // Actual insertions
    PyObject const* const range = range_location(variant.sample_start, variant.sample_end);
    if (range == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if
    return Py_BuildValue("{s:s,s:O}", "source", "observed", "location", range);
} // insertion_dict


// Build a complete variant object.
static PyObject const*
variant_dict(std::vector<mutalyzer::Variant>::const_iterator &it)
{
    PyObject const* const range = range_location(it->reference_start, it->reference_end);
    if (range == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if

    PyObject* inserted = NULL;

    // This variant starts a transposition
    if ((it->type & mutalyzer::TRANSPOSITION_OPEN) == mutalyzer::TRANSPOSITION_OPEN)
    {

        // There is always at least one inserted part of a transposition
        inserted = PyList_New(0);
        if (inserted == NULL)
        {
            Py_DECREF(range);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_New");
            return NULL;
        } // if
        PyObject const* const item = insertion_dict(*it);
        if (item == NULL)
        {
            Py_DECREF(range);
            Py_DECREF(inserted);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        if (PyList_Append(inserted, const_cast<PyObject*>(item)) != 0)
        {
            Py_DECREF(range);
            Py_DECREF(inserted);
            Py_DECREF(item);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_Append");
            return NULL;
        } // if

        // Check for a compound transposition
        if ((it->type & mutalyzer::TRANSPOSITION_CLOSE) != mutalyzer::TRANSPOSITION_CLOSE)
        {
            // As long as we do not see the end of the transposition; add these parts
            // This runs forever (until all memory is consumed) when a transposition is not closed properly.
            do
            {
                ++it;
                PyObject const* const item = insertion_dict(*it);
                if (item == NULL)
                {
                    Py_DECREF(range);
                    Py_DECREF(inserted);
                    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
                    return NULL;
                } // if
                if (PyList_Append(inserted, const_cast<PyObject*>(item)) != 0)
                {
                    Py_DECREF(range);
                    Py_DECREF(inserted);
                    Py_DECREF(item);
                    PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_Append");
                    return NULL;
                } // if
            } while ((it->type & mutalyzer::TRANSPOSITION_CLOSE) != mutalyzer::TRANSPOSITION_CLOSE);
        } // if
    } // if

    // This is not an actual variant; both reference and observed are equal.
    else if (it->type == mutalyzer::IDENTITY)
    {
        return Py_BuildValue("{s:s,s:O,s:s}", "type", "equal", "location", range, "source", "reference");
    } // if

    // This is a reverse complement variant.
    else if (it->type == mutalyzer::REVERSE_COMPLEMENT)
    {
        return Py_BuildValue("{s:s,s:O,s:s}", "type", "inversion", "location", range, "source", "reference");
    } // if

    // All other variants are deletion/insertions; guaranteed one inserted part.
    if (inserted == NULL)
    {
        inserted = PyList_New(0);
        if (inserted == NULL)
        {
            Py_DECREF(range);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_New");
            return NULL;
        } // if
        PyObject const* const item = insertion_dict(*it);
        if (item == NULL)
        {
            Py_DECREF(range);
            Py_DECREF(inserted);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        if (PyList_Append(inserted, const_cast<PyObject*>(item)) != 0)
        {
            Py_DECREF(range);
            Py_DECREF(inserted);
            Py_DECREF(item);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_Append");
            return NULL;
        } // if
    } // if
    return Py_BuildValue("{s:s,s:O,s:O,s:s}", "type", "deletion_insertion", "location", range, "inserted", inserted, "source", "reference");
} // variant_dict


// Converts the the raw C++ variants into python objects.
static PyObject*
extractor_describe_dna(PyObject*, PyObject* args)
{
    char* reference;
    char* sample;

    if (!PyArg_ParseTuple(args, "ss", &reference, &sample))
    {
        return NULL;
    } // if

    size_t const reference_length = strlen(reference);
    size_t const sample_length = strlen(sample);

    std::vector<mutalyzer::Variant> const variants = mutalyzer::extract_dna(reference, reference_length, sample, sample_length);

    PyObject* result = PyList_New(0);
    if (result == NULL)
    {
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_New");
        return NULL;
    } // if

    for (std::vector<mutalyzer::Variant>::const_iterator it = variants.begin(); it != variants.end(); ++it)
    {
        PyObject const* const item = variant_dict(it);
        if (item == NULL)
        {
            Py_DECREF(result);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        if (PyList_Append(result, const_cast<PyObject*>(item)) != 0)
        {
            Py_DECREF(item);
            Py_DECREF(result);
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for PyList_Append");
            return NULL;
        } // if
    } // for
    return result;
} // extractor_describe_dna


static PyMethodDef ExtractorMethods[] =
{
    {"describe_dna", extractor_describe_dna, METH_VARARGS,
     "Give an allele description of the change from {reference} to {observed}.\n"
     "    :arg ASCII reference: Reference sequence over the alphabet {A, C, G, T, U, (N)}.\n"
     "    :arg ASCII observed: Observed sequence over the alphabet {A, C, G, T, U, (N)}.\n"
     "    :returns list({}): A list of dictionaries representing the observed allele in terms of the reference sequence."},

    {NULL, NULL, 0, NULL}  // sentinel
}; // ExtractorMethods


static struct PyModuleDef extractormodule =
{
    PyModuleDef_HEAD_INIT,
    "description-extractor",
    "HGVS variant description extractor\n"
    "Extract a list of differences between two sequences.",
    -1,
    ExtractorMethods
}; // extractormodule


PyMODINIT_FUNC
PyInit_extractor(void)
{
    PyObject* module = PyModule_Create(&extractormodule);

    PyModule_AddStringConstant(module, "core_version", mutalyzer::VERSION);

    return module;
} // PyInit_extractor


int
main(int, char* argv[])
{
    wchar_t* const program = Py_DecodeLocale(argv[0], NULL);
    if (program == NULL)
    {
        fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
        return EXIT_FAILURE;
    } // if

    PyImport_AppendInittab("extractor", PyInit_extractor);

    Py_SetProgramName(program);

    Py_Initialize();

    PyMem_RawFree(program);

    return EXIT_SUCCESS;
} // main
