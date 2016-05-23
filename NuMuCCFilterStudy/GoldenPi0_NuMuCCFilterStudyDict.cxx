// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_NuMuCCFilterStudyDict

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "FilterEvents.h"
#include "StudyNeutrinoInteraction.h"
#include "SearchPFPartHierarchy.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLSearchPFPartHierarchy_Dictionary();
   static void larlitecLcLSearchPFPartHierarchy_TClassManip(TClass*);
   static void *new_larlitecLcLSearchPFPartHierarchy(void *p = 0);
   static void *newArray_larlitecLcLSearchPFPartHierarchy(Long_t size, void *p);
   static void delete_larlitecLcLSearchPFPartHierarchy(void *p);
   static void deleteArray_larlitecLcLSearchPFPartHierarchy(void *p);
   static void destruct_larlitecLcLSearchPFPartHierarchy(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::SearchPFPartHierarchy*)
   {
      ::larlite::SearchPFPartHierarchy *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::SearchPFPartHierarchy));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::SearchPFPartHierarchy", "SearchPFPartHierarchy.h", 25,
                  typeid(::larlite::SearchPFPartHierarchy), DefineBehavior(ptr, ptr),
                  &larlitecLcLSearchPFPartHierarchy_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::SearchPFPartHierarchy) );
      instance.SetNew(&new_larlitecLcLSearchPFPartHierarchy);
      instance.SetNewArray(&newArray_larlitecLcLSearchPFPartHierarchy);
      instance.SetDelete(&delete_larlitecLcLSearchPFPartHierarchy);
      instance.SetDeleteArray(&deleteArray_larlitecLcLSearchPFPartHierarchy);
      instance.SetDestructor(&destruct_larlitecLcLSearchPFPartHierarchy);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::SearchPFPartHierarchy*)
   {
      return GenerateInitInstanceLocal((::larlite::SearchPFPartHierarchy*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::SearchPFPartHierarchy*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLSearchPFPartHierarchy_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::SearchPFPartHierarchy*)0x0)->GetClass();
      larlitecLcLSearchPFPartHierarchy_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLSearchPFPartHierarchy_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLStudyNeutrinoInteraction_Dictionary();
   static void larlitecLcLStudyNeutrinoInteraction_TClassManip(TClass*);
   static void *new_larlitecLcLStudyNeutrinoInteraction(void *p = 0);
   static void *newArray_larlitecLcLStudyNeutrinoInteraction(Long_t size, void *p);
   static void delete_larlitecLcLStudyNeutrinoInteraction(void *p);
   static void deleteArray_larlitecLcLStudyNeutrinoInteraction(void *p);
   static void destruct_larlitecLcLStudyNeutrinoInteraction(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::StudyNeutrinoInteraction*)
   {
      ::larlite::StudyNeutrinoInteraction *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::StudyNeutrinoInteraction));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::StudyNeutrinoInteraction", "StudyNeutrinoInteraction.h", 26,
                  typeid(::larlite::StudyNeutrinoInteraction), DefineBehavior(ptr, ptr),
                  &larlitecLcLStudyNeutrinoInteraction_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::StudyNeutrinoInteraction) );
      instance.SetNew(&new_larlitecLcLStudyNeutrinoInteraction);
      instance.SetNewArray(&newArray_larlitecLcLStudyNeutrinoInteraction);
      instance.SetDelete(&delete_larlitecLcLStudyNeutrinoInteraction);
      instance.SetDeleteArray(&deleteArray_larlitecLcLStudyNeutrinoInteraction);
      instance.SetDestructor(&destruct_larlitecLcLStudyNeutrinoInteraction);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::StudyNeutrinoInteraction*)
   {
      return GenerateInitInstanceLocal((::larlite::StudyNeutrinoInteraction*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::StudyNeutrinoInteraction*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLStudyNeutrinoInteraction_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::StudyNeutrinoInteraction*)0x0)->GetClass();
      larlitecLcLStudyNeutrinoInteraction_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLStudyNeutrinoInteraction_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static TClass *larlitecLcLFilterEvents_Dictionary();
   static void larlitecLcLFilterEvents_TClassManip(TClass*);
   static void *new_larlitecLcLFilterEvents(void *p = 0);
   static void *newArray_larlitecLcLFilterEvents(Long_t size, void *p);
   static void delete_larlitecLcLFilterEvents(void *p);
   static void deleteArray_larlitecLcLFilterEvents(void *p);
   static void destruct_larlitecLcLFilterEvents(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::FilterEvents*)
   {
      ::larlite::FilterEvents *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::FilterEvents));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::FilterEvents", "FilterEvents.h", 25,
                  typeid(::larlite::FilterEvents), DefineBehavior(ptr, ptr),
                  &larlitecLcLFilterEvents_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::FilterEvents) );
      instance.SetNew(&new_larlitecLcLFilterEvents);
      instance.SetNewArray(&newArray_larlitecLcLFilterEvents);
      instance.SetDelete(&delete_larlitecLcLFilterEvents);
      instance.SetDeleteArray(&deleteArray_larlitecLcLFilterEvents);
      instance.SetDestructor(&destruct_larlitecLcLFilterEvents);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::FilterEvents*)
   {
      return GenerateInitInstanceLocal((::larlite::FilterEvents*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::FilterEvents*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLFilterEvents_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::FilterEvents*)0x0)->GetClass();
      larlitecLcLFilterEvents_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLFilterEvents_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLSearchPFPartHierarchy(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::SearchPFPartHierarchy : new ::larlite::SearchPFPartHierarchy;
   }
   static void *newArray_larlitecLcLSearchPFPartHierarchy(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::SearchPFPartHierarchy[nElements] : new ::larlite::SearchPFPartHierarchy[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLSearchPFPartHierarchy(void *p) {
      delete ((::larlite::SearchPFPartHierarchy*)p);
   }
   static void deleteArray_larlitecLcLSearchPFPartHierarchy(void *p) {
      delete [] ((::larlite::SearchPFPartHierarchy*)p);
   }
   static void destruct_larlitecLcLSearchPFPartHierarchy(void *p) {
      typedef ::larlite::SearchPFPartHierarchy current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::SearchPFPartHierarchy

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLStudyNeutrinoInteraction(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::StudyNeutrinoInteraction : new ::larlite::StudyNeutrinoInteraction;
   }
   static void *newArray_larlitecLcLStudyNeutrinoInteraction(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::StudyNeutrinoInteraction[nElements] : new ::larlite::StudyNeutrinoInteraction[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLStudyNeutrinoInteraction(void *p) {
      delete ((::larlite::StudyNeutrinoInteraction*)p);
   }
   static void deleteArray_larlitecLcLStudyNeutrinoInteraction(void *p) {
      delete [] ((::larlite::StudyNeutrinoInteraction*)p);
   }
   static void destruct_larlitecLcLStudyNeutrinoInteraction(void *p) {
      typedef ::larlite::StudyNeutrinoInteraction current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::StudyNeutrinoInteraction

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLFilterEvents(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::FilterEvents : new ::larlite::FilterEvents;
   }
   static void *newArray_larlitecLcLFilterEvents(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::FilterEvents[nElements] : new ::larlite::FilterEvents[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLFilterEvents(void *p) {
      delete ((::larlite::FilterEvents*)p);
   }
   static void deleteArray_larlitecLcLFilterEvents(void *p) {
      delete [] ((::larlite::FilterEvents*)p);
   }
   static void destruct_larlitecLcLFilterEvents(void *p) {
      typedef ::larlite::FilterEvents current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::FilterEvents

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl() {
    static const char* headers[] = {
"FilterEvents.h",
"StudyNeutrinoInteraction.h",
"SearchPFPartHierarchy.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/SOFTWARE/ROOT_v60410/include",
"/home/david/uboone/larlite/UserDev/GoldenPi0/NuMuCCFilterStudy/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$SearchPFPartHierarchy.h")))  SearchPFPartHierarchy;}
namespace larlite{class __attribute__((annotate("$clingAutoload$StudyNeutrinoInteraction.h")))  StudyNeutrinoInteraction;}
namespace larlite{class __attribute__((annotate("$clingAutoload$FilterEvents.h")))  FilterEvents;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "FilterEvents.h"
#include "StudyNeutrinoInteraction.h"
#include "SearchPFPartHierarchy.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::FilterEvents", payloadCode, "@",
"larlite::SearchPFPartHierarchy", payloadCode, "@",
"larlite::StudyNeutrinoInteraction", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_NuMuCCFilterStudy",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy() {
  TriggerDictionaryInitialization_libGoldenPi0_NuMuCCFilterStudy_Impl();
}
