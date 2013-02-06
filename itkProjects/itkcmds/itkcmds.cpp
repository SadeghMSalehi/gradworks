#include "/opt/local/include/popt.h"
#include "itkinfo.h"
#include "stdio.h"
#include "stdlib.h"

extern "C" {
    /* Data values for the options. */
    static int intVal = 55;
    static int print = 0;
    static char* stringVal;
    void callback(poptContext context,
                  enum poptCallbackReason reason,
                  const struct poptOption * option,
                  const char * arg,
                  const void * data)
    {
        switch(reason)
        {
            case POPT_CALLBACK_REASON_PRE:
                printf("\t Callback in pre setting\n"); break;
            case POPT_CALLBACK_REASON_POST:
                printf("\t Callback in post setting\n"); break;
            case POPT_CALLBACK_REASON_OPTION:
                printf("\t Callback in option setting\n"); break;
        }
    }

    static struct poptOption optionsTable[] = {
        { "int", (char) 'i', POPT_ARG_INT, (void*) &intVal, 0,
            "follow with an integer value", "2, 4, 8, or 16" },
        {  "file", (char) 'f', POPT_ARG_STRING, (void*) &stringVal, 0,
             "follow with a file name", NULL },
        {  "print", (char) 'p', POPT_ARG_NONE, &print, 0,
             "send output to the printer", NULL },
        POPT_AUTOALIAS
        POPT_AUTOHELP
        POPT_TABLEEND
    };
}

int main(int argc, const char* argv[]) {
    mainImageInfo(argc, argv);
}

int main_poptTest(int argc, const char *argv[]) {
    poptContext context = poptGetContext("itkcmds", argc, argv, optionsTable, 0);
    int option = -1;

    while ((option = poptGetNextOpt(context)) > 0) {
        printf("option = %d\n", option);
    };
    const char** args = poptGetArgs(context);
    argc = 0;
    while (args[argc++] != NULL) {
        printf("%s\n", args[argc-1]);
    }

    poptPrintUsage(context, stdout, 2);

    printf("option = %d\n", option);
    /* Print out option values. */
    printf("After processing, options have values:\n");
    printf("\t intVal holds %d\n", intVal);
    printf("\t print flag holds %d\n", print);
    printf("\t stringVal holds [%s]\n", stringVal);
    poptFreeContext(context);
    printf("Remaining arguments: %d\n", argc);
    exit(0);
}