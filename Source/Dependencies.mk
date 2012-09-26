#
#  Dependencies:
#
$(LIBRARY)(ConstantsModule.o):    $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(FilesModule.o):        $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(StringsModule.o):      $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(StringsModule.o):      $(LIBRARY)(ConstantsModule.o)
$(LIBRARY)(CommandModule.o):      $(LIBRARY)(StringsModule.o)
$(LIBRARY)(CommandModule.o):      $(LIBRARY)(FilesModule.o)
$(LIBRARY)(OptimModule.o):        $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(OptimModule.o):        $(LIBRARY)(ConstantsModule.o)
#
