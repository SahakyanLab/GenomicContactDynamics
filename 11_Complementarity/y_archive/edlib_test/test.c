#include <stdio.h>
#include "edlib.h"

int main() {
    EdlibAlignResult result = edlibAlign("TELEPHONE", 9, "ALEAAONPA", 9, edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status == EDLIB_STATUS_OK) {
        printf("%d\n", result.editDistance);
        printf("%d\n", result.alignmentLength);
        printf("%d\n", result.endLocations[0]);
    }
    
    char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
    printf("%s", cigar);
    
    edlibFreeAlignResult(result);
}
