use  bindgen;
use cc;

use std::{env, path::PathBuf};

fn main() {
    cc::Build::new().file("lib/hello.c").compile("hello");
    let bindings = bindgen::Builder::default()
        .header("lib/hello.h")
        .generate()
        .expect("Unable to generate bindings");

    let out_dir = PathBuf::from(env::var("OUT_DIR").unwrap());
    bindings
        .write_to_file(out_dir.join("bindings.rs"))
        .expect("Couldn't write bindings!");
}
