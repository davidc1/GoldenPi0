// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME GoldenPi0_MakeMCVertexDict

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
#include "MakeMCVertex.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static TClass *larlitecLcLMakeMCVertex_Dictionary();
   static void larlitecLcLMakeMCVertex_TClassManip(TClass*);
   static void *new_larlitecLcLMakeMCVertex(void *p = 0);
   static void *newArray_larlitecLcLMakeMCVertex(Long_t size, void *p);
   static void delete_larlitecLcLMakeMCVertex(void *p);
   static void deleteArray_larlitecLcLMakeMCVertex(void *p);
   static void destruct_larlitecLcLMakeMCVertex(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::larlite::MakeMCVertex*)
   {
      ::larlite::MakeMCVertex *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::larlite::MakeMCVertex));
      static ::ROOT::TGenericClassInfo 
         instance("larlite::MakeMCVertex", "MakeMCVertex.h", 26,
                  typeid(::larlite::MakeMCVertex), DefineBehavior(ptr, ptr),
                  &larlitecLcLMakeMCVertex_Dictionary, isa_proxy, 4,
                  sizeof(::larlite::MakeMCVertex) );
      instance.SetNew(&new_larlitecLcLMakeMCVertex);
      instance.SetNewArray(&newArray_larlitecLcLMakeMCVertex);
      instance.SetDelete(&delete_larlitecLcLMakeMCVertex);
      instance.SetDeleteArray(&deleteArray_larlitecLcLMakeMCVertex);
      instance.SetDestructor(&destruct_larlitecLcLMakeMCVertex);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::larlite::MakeMCVertex*)
   {
      return GenerateInitInstanceLocal((::larlite::MakeMCVertex*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::larlite::MakeMCVertex*)0x0); R__UseDummy(_R__UNIQUE_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *larlitecLcLMakeMCVertex_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const ::larlite::MakeMCVertex*)0x0)->GetClass();
      larlitecLcLMakeMCVertex_TClassManip(theClass);
   return theClass;
   }

   static void larlitecLcLMakeMCVertex_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_larlitecLcLMakeMCVertex(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::MakeMCVertex : new ::larlite::MakeMCVertex;
   }
   static void *newArray_larlitecLcLMakeMCVertex(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::larlite::MakeMCVertex[nElements] : new ::larlite::MakeMCVertex[nElements];
   }
   // Wrapper around operator delete
   static void delete_larlitecLcLMakeMCVertex(void *p) {
      delete ((::larlite::MakeMCVertex*)p);
   }
   static void deleteArray_larlitecLcLMakeMCVertex(void *p) {
      delete [] ((::larlite::MakeMCVertex*)p);
   }
   static void destruct_larlitecLcLMakeMCVertex(void *p) {
      typedef ::larlite::MakeMCVertex current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::larlite::MakeMCVertex

namespace {
  void TriggerDictionaryInitialization_libGoldenPi0_MakeMCVertex_Impl() {
    static const char* headers[] = {
"MakeMCVertex.h",
0
    };
    static const char* includePaths[] = {
"/home/david/uboone/larlite/core",
"/home/david/SOFTWARE/ROOT_v60410/include",
"/home/david/uboone/larlite/UserDev/GoldenPi0/MakeMCVertex/",
0
    };
    static const char* fwdDeclCode = 
R"DICTFWDDCLS(
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace larlite{class __attribute__((annotate("$clingAutoload$MakeMCVertex.h")))  MakeMCVertex;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "MakeMCVertex.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"larlite::MakeMCVertex", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("libGoldenPi0_MakeMCVertex",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_libGoldenPi0_MakeMCVertex_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_libGoldenPi0_MakeMCVertex_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_libGoldenPi0_MakeMCVertex() {
  TriggerDictionaryInitialization_libGoldenPi0_MakeMCVertex_Impl();
}
