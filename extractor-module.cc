#include <Python.h>

#include "extractor-core/include/extractor.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

#include <vector>


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
        PyObject* item = Py_BuildValue("{s:i, s:i, s:i, s:i, s:i, s:i, s:i}",
                                       "reference_start", it->reference_start,
                                       "reference_end", it->reference_end,
                                       "sample_start", it->sample_start,
                                       "sample_end", it->sample_end,
                                       "type", it->type,
                                       "transposition_start", it->transposition_start,
                                       "transposition_end", it->transposition_end);
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
