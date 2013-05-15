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
$(LIBRARY)(LibF95.o):  $(LIBRARY)(IntrinsicTypesModule.o)
$(LIBRARY)(LibF95.o):  $(LIBRARY)(ConstantsModule.o)
$(LIBRARY)(LibF95.o):  $(LIBRARY)(FilesModule.o)
$(LIBRARY)(LibF95.o):  $(LIBRARY)(OptimModule.o)
$(LIBRARY)(LibF95.o):  $(LIBRARY)(TimerModule.o)
$(LIBRARY)(LibF95.o):  $(LIBRARY)(StringsModule.o)
$(LIBRARY)(LibF95.o):  $(LIBRARY)(CommandModule.o)
