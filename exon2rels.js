#!/usr/bin/env node
var fs     = require("fs");
var cl     = require("termcolor").define();
var LS     = require("linestream");
var Junjo  = require("junjo");


function main() {
  var $j = new Junjo({noTimeout : true});

  $j("filename", function(filename) {
    console.assert(require("path").existsSync(filename));
    return filename;
  });

  $j("exons", function(filename) {
    return fs.readFile(filename, "utf-8", this.cb);
  })
  .post(function(data) {
    return data.split("\n");
  })
  .after("filename")
  .eshift();

  $j("groups", function() {
    process.stdin.setEncoding("utf8");
    process.stdin.resume();

    var lines = new LS(process.stdin, {trim: true});
    this.absorb(lines, "data", function(line, ret) {
      if (!ret) ret = {};
      ret[line] = true;
      return ret;
    });
  })
  .eshift();

  $j("check", function(exons, groups) {
    this.$.exons = exons;
    this.iterate(exons, function(ex, k) {
      var data = ex.split("\t");
      if (data.length < 10) return;

      var startCount = Number(data[5]);
      if (!startCount) return;

      var strand = data[3];

      var starts = data[9].split(",").filter(function(v) {
        return v.match(/\[1\//);
      })
      .map(function(v) {
        return v.split("[")[0];
      });

      var delta = (strand == "+") ? 1 : -1;

      starts.forEach(function(name, k2) {
        var  current = k - delta;
        while (true) {
          current += delta;
          var ex = exons[current];
          if (!ex) return;
          var data = ex.split("\t");
          if (data.length < 11) continue;
          var thisStrand = data[3];
          if (strand != thisStrand) continue;

          var splices = data[9].split(",");

          var hasTerm = splices.some(function(v) {
            var result = v.match(/\[([0-9]+)\/([0-9]+)\]/);
            var thisname = v.split("[")[0];
            return result && result[1] == result[2] && thisname == name;
          });
          if (hasTerm) { return }

          var splicenames = splices.map(function(v) {
            return v.split("[")[0];
          });

          if (!splicenames.length) continue;

          splicenames.forEach(function(thisname, i) {
            if (thisname != name) {
              var keyarr = [thisname, name].sort();
              if (!groups[keyarr.join(",")]) {
                console.log(keyarr.join("\t"));
              }
            }
          }, this);
        }
      }, this);
    });
  })
  .after("exons", "groups")

  $j.run(process.argv[2]);
}


if (process.argv[1].match('/([^/]+?)(\.js)?$')[1] == __filename.match('/([^/]+?)(\.js)?$')[1]) main();
