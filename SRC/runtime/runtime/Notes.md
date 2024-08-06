

- `LibraryBroker<type>`


```
  static library = G3_Library(
     LibraryBroker<UniaxialMaterial>(&uniaxial_material_library),
     LibraryBroker<NDMaterial>(&multiaxial_material_library),
  );

read
parse
build
define
broker
invoke   call
create   make             new
insert   save             add
delete   wipe             del
update    -
remove   wipe   clean
         read   build     get
         scan 
         show
         load
         dump
                clear
                count
                erase
                print


del
add
get
set
put
log
let
new

->callObj(tag, task, resp) -> state
->callObj(obj, task, resp)
->makeObj(tag)
->readObj(argc, argv) -> Obj*

->saveObj(obj)
->wipeObj(tag)

->pullObj(obj)
->yankObj(obj)


  static broker = G3_Broker(library);

  T* library->readNew<X>(builder*, argc, argv);
  library->makeNew<X>();


newObj(tag)
getObj(argc, argv);
addObj(getObj(argc, argv));


  G3_RuntimeBroker(library)




  G3_ModelBuilder()
  model_builder->makeNew<X>(inst_tag)
  model_builder->tag<X>(X*)
  model_builder->readNew(argc, argv);


  runtime_broker->newEmpty<X>(class_tag);
```




ElemRoutine(rt, state, task, argc, argv)








