/*-----------------------------------------------------------------*/
/* Local functions for fragmentation interface                     */
/*-----------------------------------------------------------------*/
/* copy-input.c                                                    */
void copySimParam(GENERAL_DATA *,BONDED *,CLASS *,CP *,GENERAL_DATA *,BONDED *,
                  CLASS *,CP *,CLASS_PARSE *,CP_PARSE *,FILENAME_PARSE *);

/*-----------------------------------------------------------------*/
/* control-mol-params-frag.c                                       */

void controlMolParamsFrag(CLASS *,GENERAL_DATA *,BONDED *,CP *,CLASS *,
                          GENERAL_DATA *,BONDED *,CP *,CLASS_PARSE *,
                          CP_PARSE *,FREE_PARSE *,FILENAME_PARSE *);

void controlSetMolParamsFrag(CP_PARSE *,CLASS_PARSE *,
                FILENAME_PARSE *,FREE_PARSE *,BONDED *,CLASS *,CP *,
                GENERAL_DATA *,CLASS *class,CP *,GENERAL_DATA *,DICT_MOL *,DICT_WORD *,
                char *,int *,int ,double ,int ,int ,DICT_MOL *);



