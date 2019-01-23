#include <Python.h>

#include "extractor-core/include/extractor.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>


static PyObject*
point_location(size_t const point)
{
    return Py_BuildValue("{s:s,s:i}", "type", "point", "position", point);
} // point_location


static PyObject*
range_location(size_t const start, size_t const end)
{
    PyObject* start_object = point_location(start);
    PyObject* end_object = point_location(end);

    if (start_object == NULL || end_object == NULL)
    {
        Py_DECREF(start_object);
        Py_DECREF(end_object);
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if
    return Py_BuildValue("{s:s,s:O,s:O}", "type", "range", "start", start_object, "end", end_object);
} // range_location


static PyObject*
variant_dict(mutalyzer::Variant const &variant)
{
    PyObject* reference_range = range_location(variant.reference_start, variant.reference_end);
    PyObject* sample_range = range_location(variant.sample_start, variant.sample_end);

    if (reference_range == NULL || sample_range == NULL)
    {
        Py_DECREF(reference_range);
        Py_DECREF(sample_range);
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if

    if (variant.type == mutalyzer::IDENTITY)
    {
        return Py_BuildValue("{s:s,s:O}", "type", "equal", "location", sample_range);
    } // if
    else if (variant.type == mutalyzer::REVERSE_COMPLEMENT)
    {
        return Py_BuildValue("{s:s,s:O}", "type", "inv", "location", sample_range);
    } // if

    PyObject* inserted = Py_BuildValue("{s:s,s:O}", "source", "observed", "location", reference_range);
    if (inserted == NULL)
    {
        Py_DECREF(reference_range);
        Py_DECREF(sample_range);
        PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
        return NULL;
    } // if

    return Py_BuildValue("{s:s,s:O,s:[O]}", "type", "delins", "location", sample_range, "insertions", inserted);
} // variant_dict


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
        // TODO: transpositions
        PyObject* item = variant_dict(*it);
        if (item == NULL)
        {
            PyErr_SetString(PyExc_MemoryError, "Could not allocate memory for Py_BuildValue");
            return NULL;
        } // if
        if (PyList_Append(result, item) != 0)
        {
            Py_DECREF(item);
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
    "   :arg ascii reference: Reference sequence.\n"
    "   :arg ascii observed: Observed sequence.\n"
    "   :returns list({}): A list of dictionaries representing the obsereved allele in terms of the reference sequence."
    },

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
