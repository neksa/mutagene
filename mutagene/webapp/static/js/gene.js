
function resetScanSeqResults() {
  $("#gene_results_labels_placeholder").empty();
  $("#gene_results_placeholder").empty();
}


function drawScanScatterPlot(json) {
    if (json.error) {
      console.log('AJAX request returned error: ' + json.error)
      return;
    }
    
    console.log(json);

    $("#scanResultsModal").modal('show');

    ref = "#scan_results_scatterplot";
    $("#scan_results_scatterplot_loading").text("");
    d3.selectAll(ref + " > *").remove();
    var data = json.data;

    var dot_size = 4.5;
    var width = 500;
    var height = 430;
    var margin = 80;
    var svg = d3.select(ref)
        .attr("width", width + 2 * margin)
        .attr("height", height + margin)
      .append("g")
        .attr("transform", "translate(" + (margin) + "," + 10 + ")");

    var x = d3.scale.log().domain([1, d3.max(data, function(d){return d.observed;})]).nice().range([20, width]);
    var y = d3.scale.log().domain([d3.min(data, function(d){ return d.mutability;}), d3.max(data, function(d){return d.mutability;})]).nice().range([height, 0]);
    var xAxis = d3.svg.axis().scale(x).orient("bottom").ticks(5);  // .tickFormat(d3.format('d'));
    var yAxis = d3.svg.axis().scale(y).orient("left").ticks(5);  // .tickFormat(d3.format('d'));

    var color = d3.scale.linear().domain([0.0, 0.5, 1.0]).range(['#1a9641', "#ffffbf", '#d7191c']);
    console.log(color(0.0));
    console.log(color(0.5));
    console.log(color(1.0));

    // var min_val = Math.min(d3.min(data, yValue), d3.min(data, xValue)) - 0.01;
    // var max_val = Math.max(d3.max(data, yValue), d3.max(data, xValue)) + 0.01;
    // yScale.domain([0, d3.max(data, function(d){d.})]);
    // xScale.domain([0, d3.max(data, xValue)]);

    // x-axis
    svg.append("g")
        .attr("class", "x axis")
        .attr("transform", "translate(0," + height + ")")
        .call(xAxis)
      .append("text")
        .attr("class", "scatterplot-label")
        .attr("x", width)
        .attr("y", 35)
        .style("text-anchor", "end")
        .text('Frequency of observed cancer mutations (log scale)');

    // y-axis
    svg.append("g")
        .attr("class", "y axis")
        .call(yAxis)
      .append("text")
        .attr("class", "scatterplot-label")
        .attr("transform", "rotate(-90)")
        .attr("y", -55)
        .attr("dy", ".71em")
        // .attr("dy", "1em")
        .style("text-anchor", "end")
        .text('Expected mutability according to mutational profile or signature (log scale)');

    var tooltip = d3.select(".scatterplot-tooltip");

    svg.selectAll(".dot")
         .data(data)
           .enter().append("circle")
             .attr('class', 'scatter-dot')
             .attr('r', dot_size)
             .attr('opacity', 0.7)
             .attr('cx', function(d){ return x(d.observed + Math.random()*0.25 - 0.25/2)})
             .attr('cy', function(d){ return y(d.mutability); })
             .style('fill', function(d) {
                  var impact = parseFloat(d.impact);
                  if (impact > 0) {
                      return color(impact);
                  } else {
                      return "yellow";
                  }
             })
             .on("mouseover", function(d) {
                 tooltip.transition()
                      .duration(300)
                      .style("opacity", .95);
                      // .style('width', 300);
                 impact = "unknown";
                 if (d.impact > 0) {
                   impact = d.impact.toFixed(2);
                 }
                 tooltip.html(d.mutation + "&nbsp;(FATHMM&nbsp;score&nbsp;=&nbsp;" + impact + ")")
                      .style("left", (d3.event.pageX + 5) + "px")
                      .style("top", (d3.event.pageY - 28) + "px");
             })
             .on("mouseout", function(d) {
                 tooltip.transition()
                      .duration(500)
                      .style("opacity", 0);
             });

    // var show_legend = true;
    // if (show_legend) {
    
    //   var legend = svg.selectAll(".legend")
    //       .data(mutations.slice())
    //     .enter().append("g")
    //       .attr("class", "legend")
    //       .attr("transform", function(d, i) { return "translate(0," + (i * 16 + 10) + ")"; });

    //   legend.append("circle")
    //       .attr("cx", width+margin - 10)
    //       .attr('r', 6)
    //       .attr('class', 'scatter-dot')
    //       .style("fill", color);

    //   legend.append("text")
    //       .attr("x", width+margin - 20)
    //       .attr("y", 0)
    //       .attr("dy", 4)
    //       .attr("font-family", "Courier")
    //       .attr("font-size", "12px")
    //       .style("text-anchor", "end")
    //       .text(function(d) { return d[0] + ' → ' + d[1]; });

    //   svg.append("text")
    //       .attr("x", width + margin - 0)
    //       .attr("y", 0)
    //       .style("text-anchor", "end")
    //       .text("Mutation:");
    // }
}

var firstTime = true;

function drawScanSeqResults(json) {
    if (json.error) {
      console.log('AJAX request returned error: ' + json.error)
      return;
    }
    // console.log(json);
    // console.log(json.na_records);
    // console.log(json.aa_records);
    // console.log(json.aa_matrix);
    $(".results_legend").show();
    if(firstTime) {
        firstTime = false;
        scrollToAnchor("results");
    }

    var x_offset        = 10;

    var y_offset_observed = 0;  // baseline
    var y_height_observed = 90;
    var y_offset_dna      = 0; // baseline
    var y_height_dna      = 100;
    var y_height_na_matrix = 75;
    var y_offset_aa       = 100 + y_height_na_matrix;  // baseline
    var y_height_aa       = 75;
    var y_offset_matrix   = 210 + y_height_na_matrix; // baseline
    var y_offset_na_matrix= y_offset_dna + y_height_dna + 15; // baseline

    // var dna_width = 15;
    var dna_width = 8;
    var spacer = 1;
    var aa_width = 3 * dna_width + 2 * spacer;
    // var matrix_cell_height = 15;
    var matrix_cell_height = 14;

    dna_X = json.na_records.length;
    dna_Y = json.nucleic_acids.length;
    aa_X = json.aa_records.length;
    aa_Y = json.amino_acids.length;

    var y_height_na_matrix = dna_Y * (matrix_cell_height + spacer);
    var y_height_matrix = aa_Y * (matrix_cell_height + spacer);

    var plot_height = y_offset_matrix + y_height_matrix;
    var plot_width = x_offset + dna_X * (dna_width + spacer);

    var labels_plot_height = plot_height;
    var labels_plot_width = 150; // 116;
    var labels_right_border = labels_plot_width - 49// - 15;
    // var data = json.aa_matrix;
    // var matrix_dim = {
    //   X: json.aa_records.length,
    //   Y: 20
    // };
    // var matrix_cell = {w: 11, h: 11, spacer: 1};

    // var dna_offset = {
    //   x:,
    //   y:
    // };

    // var aa_offset = {
    //   x: dna_offset.x,
    //   y: dna_offset.y + 100
    // };

    // var matrix_offset = {
    //     x: aa_offset.x,
    //     y: aa_offset.y + 100
    //   };

    // var aa_mutabilities = function(d) { return d.x;}
    // var na_mutabilities = function(d) { return d.x;}
    var na_max = d3.max(json.na_records, function(d) { return parseFloat(d.mutability); });
    var aa_max = d3.max(json.aa_records, function(d) { return parseFloat(d.mutability); });
    var na_matrix_max = d3.max(json.na_matrix, function(d) { return parseFloat(d.freq); });
    var aa_matrix_max = d3.max(json.aa_matrix, function(d) { return parseFloat(d.freq); });

    var na_line_scale = d3.scale.linear()
        .domain([0, na_max])
        .range([y_height_dna - 10, 3]);

    // console.log("na_max", na_max);
    // console.log("aa_max", aa_max);
    // console.log("aa_matrix_max", aa_matrix_max);
    // var aa_matrix_scale = d3.scale.linear().domain([0, aa_matrix_max]).range([0, 1.0]);
    //var colorScaleDNA = d3.scale.linear().domain([0, na_matrix_max]).range(["#F0F0F0",  "#b011b0"]);

    // var colorScaleMissense = d3.scale.linear().domain([0, aa_matrix_max]).range(["#F0F0F0",  "#713b85"]);
    var colorScaleMissense = d3.scale.linear().domain([0, aa_matrix_max]).range(["#F4F0F4",  "#713b85"]);
    // var colorScaleSilent = d3.scale.linear().domain([0, aa_matrix_max]).range(["#F0F0F0",  "#50853b"]);
    var colorScaleSilent = d3.scale.linear().domain([0, aa_matrix_max]).range(["#EBFCE5",  "#50853b"]);
    // var colorScaleNonsense = d3.scale.linear().domain([0, aa_matrix_max]).range(["#F0F0F0",  "#bc2819"]);
    var colorScaleNonsense = d3.scale.linear().domain([0, aa_matrix_max]).range(["#F6DEDC",  "#bc2819"]);

    resetScanSeqResults();

    var labels_id = "#gene_results_labels_placeholder";
    // $(labels_id).empty();
    var svg_labels = d3.select(labels_id)
        .append("svg")
        .attr("width", labels_plot_width)
        .attr("height", labels_plot_height);

    svg_labels.append("rect")
        .attr("x", "0")
        .attr("y", "0")
        .attr("width", labels_plot_width)
        .attr("height", labels_plot_height)
        .attr("fill", "rgba(240, 240, 240, 0.5)");

    var placeholder_id = '#gene_results_placeholder';
    // $(placeholder_id).empty();
    // d3.selectAll("#vis > *").remove();
    var svg = d3.select(placeholder_id)
        .append("svg")
        .attr("width", plot_width)
        .attr("height", plot_height);

    // init tooltip
    var tip = d3.tip()
      .attr('class', 'd3-tip')
      .offset([-10, 0])
      .html(function(d) {
          var gene = $("#gene").val();
          var html = ""
          // alert(d.seq_num);
          var mutation = "mutations";
          if (d.observed_num == 1) {
              mutation = "mutation";
          }
          // TIP works for na_records and aa_matrix. na_records have seq_num
          if (d.seq_num === undefined) {
              if (d.aa_from_long === undefined) {
                  html += "DNA substitution " + d.na_from + " → " + d.na_to + "<br>Mutability = " + d3.format(".2f")(d.freq) + " per Megabase<br>"; //  .toExponential(2);
                
              } else {
                  html += "Amino acid substitution " + d.aa_from_long + " → " + d.aa_to_long + "<br>Mutability = " + d3.format(".2f")(d.freq) + " per Megabase<br>"; //  .toExponential(2);
              }
              if (d.observed_num > 0) {
                  html += "<br>Number of observed mutations in:"
                  for (k in d.observed) {
                      html += "<br>&nbsp;&nbsp;&nbsp;&nbsp;" + k + ": " + d.observed[k];
                  }
              }
          } else {
              if (d.observed_num > 0) {
                  html += "Number of observed mutations in " + gene + " on DNA position " + (d.seq_num + 1) + ":<br>";
                  for (k in d.observed) {
                      html += "<br> " + k + ": " + d.observed[k];
                  }
                  // html += "Observed " + d.observed + " pan-cancer " + mutation + " in DNA position " + (d.seq_num + 1);
              }
          }
          return html;
      });
    svg.call(tip);

    //////////// baselines ///////
    svg.append("line")
          .attr("x1", 0)
          .attr("y1", 0)
          .attr("x2", plot_width)
          .attr("y2", 0)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", 0)
          .attr("x2", labels_plot_width)
          .attr("y2", 0)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    svg.append("line")
          .attr("x1", 0)
          .attr("y1", plot_height)
          .attr("x2", plot_width)
          .attr("y2", plot_height)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", labels_plot_height)
          .attr("x2", labels_plot_width)
          .attr("y2", labels_plot_height)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    svg.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_observed)
          .attr("x2", plot_width)
          .attr("y2", y_offset_observed)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_observed)
          .attr("x2", labels_plot_width)
          .attr("y2", y_offset_observed)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    svg.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_dna)
          .attr("x2", plot_width)
          .attr("y2", y_offset_dna)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_dna)
          .attr("x2", labels_plot_width)
          .attr("y2", y_offset_dna)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    svg.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_aa)
          .attr("x2", plot_width)
          .attr("y2", y_offset_aa)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_aa)
          .attr("x2", labels_plot_width)
          .attr("y2", y_offset_aa)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    svg.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_matrix)
          .attr("x2", plot_width)
          .attr("y2", y_offset_matrix)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_matrix)
          .attr("x2", labels_plot_width)
          .attr("y2", y_offset_matrix)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    svg.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_na_matrix)
          .attr("x2", plot_width)
          .attr("y2", y_offset_na_matrix)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");
    svg_labels.append("line")
          .attr("x1", 0)
          .attr("y1", y_offset_na_matrix)
          .attr("x2", labels_plot_width)
          .attr("y2", y_offset_na_matrix)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2");

    // vertical line
    svg_labels.append("line")
          .attr("x1", labels_right_border + 7)
          .attr("y1", 0)
          .attr("x2", labels_right_border + 7)
          .attr("y2", labels_plot_height)
          .attr("stroke", "rgba(150, 150, 150, 0.5)")
          .attr("stroke-width", "2")


    svg.selectAll(".grid")
        .data(json.aa_records)
        .enter().append('line')
        // .attr('y', function(d, i) {
        //     return y_offset_aa;
        // })
        .attr('x1', function(d, i) {
            var x = d.seq_num + 1;
            return x_offset + x * (aa_width + spacer);
        })
        .attr('x2', function(d, i) {
            var x = d.seq_num + 1;
            return x_offset + x * (aa_width + spacer);
        })
        .attr('class', 'grid-line')
        .attr("y1", 0)
        .attr("y2", plot_height);
        // .attr("stroke", "rgba(200, 200, 200, 0.5)")
        // .attr("stroke-width", "1");

    /////////////////////////////////////////////////////////  AXES //////////////////////////////////////////////////////////////////
    var na_axis = d3.svg.axis()
        .scale(na_line_scale)
        .orient("left")
        .ticks(3);
        // .tickFormat(d3.format("e"));
    
    svg_labels.append("g")
        .attr("class", "mutability-axis")
        .attr("transform", "translate(" + (labels_plot_width - 2) + ",0)")
        .style("font-size", '10px')
        .style("font-weight", 'normal')
        .style("font-family", 'sans-serif')
        .call(na_axis);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    var plotCell = svg.selectAll(".matrix")
         // .selectAll(".cell")
         .data(json.na_matrix.filter(function(d){
              var freq = +d.freq || 0;
              return freq > 0.0;
           }))
           .enter()
           .append("svg:rect")
             // .attr('class', 'gcell')
             .attr('width', dna_width)
             .attr('height', matrix_cell_height)
             .attr('x', function(d) {
                return d.x * (dna_width + spacer) + x_offset;
             })
             .attr('y', function(d) {
                return y_offset_na_matrix + d.y * (matrix_cell_height + spacer);
             })
             .style('fill', function(d) {
                var freq = +d.freq || 0;
                if (d.result == 'S') {
                  return colorScaleSilent(freq);
                }
                
                if (d.result == 'M') {
                  return colorScaleMissense(freq);
                }
                
                if (d.result == 'N') {
                  return colorScaleNonsense(freq);
                }
             })
             .on('click', function(d) {
                // console.log('clicked on a cell ', d);
                var text = "";
                for (k in d.observed) {
                   text += k + ": " + d.observed[k] + "\n";
                }
                alert(text);
                // alert('clicked on a cell ', d.pos);
             })
             .on('mouseover', tip.show)
             .on('mouseout', tip.hide); 

    var plotCell = svg.selectAll(".matrix")
         // .selectAll(".cell")
         .data(json.aa_matrix.filter(function(d){
              var freq = +d.freq || 0;
              return freq > 0.0;
           }))
           .enter()
           .append("svg:rect")
             // .attr('class', 'gcell')
             .attr('width', aa_width)
             .attr('height', matrix_cell_height)
             .attr('x', function(d) {
                return d.x * (aa_width + spacer) + x_offset;
             })
             .attr('y', function(d) {
                return y_offset_matrix + d.y * (matrix_cell_height + spacer);
             })
             .style('fill', function(d) {
                var freq = +d.freq || 0;
                var P = d.aa_from;
                var Q = d.aa_to;
                if (P == Q) {
                  return colorScaleSilent(freq);
                }
                if ((P != "*" && Q == "*") || (P == "*" && Q != "*") || (d.pos == 1 && P != Q)) {
                  return colorScaleNonsense(freq);
                }
                return colorScaleMissense(freq);
             })
             .on('click', function(d) {
                // console.log('clicked on a cell ', d);
                var text = "";
                for (k in d.observed) {
                   text += k + ": " + d.observed[k] + "\n";
                }
                alert(text);
                // alert('clicked on a cell ', d.pos);
             })
             .on('mouseover', tip.show)
             .on('mouseout', tip.hide); 


    /////////////////////////////////  DNA SEQUENCE NUMBERS ///////////////////
    var dnaSequenceNum = svg.append("g")
        .attr("class", "dna-indices")
        .selectAll(".dna-indices")
        
        // .data(json.na_records)
        
        .data(json.na_records.filter(function(d, i) { return (i + 1) % 3 == 1; }))

        .enter().append('text')
        .attr('y', function(d, i) {
            return y_offset_dna + y_height_dna + 12;
        })
        .attr('x', function(d, i) {
            return x_offset + d.seq_num * (dna_width + spacer) + 2;
        })
        .text(function(d, i) {
            return i*3 + 1;
        });


    /////////////////////////////////  AMINO ACID SEQUENCE TEXT ///////////////////
    var aaSequence = svg.append("g")
        .attr("class", "aa-sequence")
        .selectAll(".aa-sequence")
        .data(json.aa_records)
        .enter().append('text')
        .attr('y', function(d, i) {
            // return y_offset_aa + y_height_aa;
            return y_offset_aa + y_height_aa + 30;
        })
        .attr('x', function(d, i) {
            return x_offset + d.seq_num * (aa_width + spacer) + aa_width / 2;
        })
        .text(function(d, i) {
            return d.seq_long;
        });


    /////////////////////////////////  AMINO ACID SEQUENCE NUMBERS ///////////////////
    var aaSequenceNum = svg.append("g")
        .attr("class", "aa-indices")
        .selectAll(".aa-indices")
        .data(json.aa_records)
        .enter().append('text')
        .attr('y', function(d, i) {
            // return y_offset_aa + 15;
            return y_offset_aa + y_height_aa + 20;
        })
        .attr('x', function(d, i) {
            return x_offset + d.seq_num * (aa_width + spacer) + aa_width / 2;
        })
        .attr('title', function(d) {return d.seq_num + 1;})
        .text(function(d, i) {
            return d.seq_num + 1;
        });


    var naLabels = svg_labels.append("g")
        .attr("class", "rowLabels")
        .selectAll(".rowLabels")
        .data(json.nucleic_acids)
        .enter().append('text')
        .attr('y', function(d, i) {
            // console.log("y",  y_offset_matrix + (i + 1) * (matrix_cell_height + spacer))
            return y_offset_na_matrix + (i) * (matrix_cell_height + spacer) + matrix_cell_height * 0.8;
        })
        .attr('x', function(d, i) {
            return labels_right_border;
        })
        .attr('title', function(d) {return d;})
        .text(function(d, i) {
            return d;
        });

    var aaLabels = svg_labels.append("g")
        .attr("class", "rowLabels")
        .selectAll(".rowLabels")
        .data(json.amino_acids)
        .enter().append('text')
        .attr('y', function(d, i) {
            // console.log("y",  y_offset_matrix + (i + 1) * (matrix_cell_height + spacer))
            return y_offset_matrix + (i) * (matrix_cell_height + spacer) + matrix_cell_height * 0.8;
        })
        .attr('x', function(d, i) {
            return labels_right_border;
        })
        .attr('title', function(d) {return d;})
        .text(function(d, i) {
            return d;
        });

    svg_labels.append("text")
        .attr("x", labels_right_border)
        .attr("y", y_offset_aa + y_height_aa / 4)
        .attr("dy", "1.4em")
        .style("text-anchor", "end")
        .text("Amino acid");
    svg_labels.append("text")
        .attr("x", labels_right_border)
        .attr("y", y_offset_aa + y_height_aa / 4)
        .attr("dy", "2.8em")
        .style("text-anchor", "end")
        .text("mutability");
    svg_labels.append("text")
        .attr("x", labels_right_border)
        .attr("y", y_offset_aa + y_height_aa / 4)
        .attr("dy", "4.2em")
        .style("text-anchor", "end")
        .text("(per Megabase)");

    svg_labels.append("text")
        .attr("x", labels_right_border)
        .attr("y", y_offset_dna + y_height_dna / 4)
        .attr("dy", "1.4em")
        .style("text-anchor", "end")
        .text("DNA");
    svg_labels.append("text")
        .attr("x", labels_right_border)
        .attr("y", y_offset_dna + y_height_dna / 4)
        .attr("dy", "2.8em")
        .style("text-anchor", "end")
        .text("mutability");
    svg_labels.append("text")
        .attr("x", labels_right_border)
        .attr("y", y_offset_dna + y_height_dna / 4)
        .attr("dy", "4.2em")
        .style("text-anchor", "end")
        .text("(per Megabase)");

    // //////////////////////////  OBSERVED DNA MUTATIONS BARS //////////////////////////////////////////////////////////////
    // svg.selectAll("bars")
    //           .data(json.na_records)
    //       .enter().append("svg:rect")
    //       .filter(function(d){ return d.observed_num > 0;})
    //           // .attr("class", "line")
    //           .attr("stroke", "black")
    //           .attr("stroke-width", 0.5)
    //           .attr("height", function(d) {
    //                 return (Math.round(Math.log(d.observed_num) + 1)*14);
    //              })
    //           .attr("width", dna_width)
    //           .attr("x", function(d) { return x_offset + d.seq_num * (dna_width + spacer); })
    //           .attr("y", function(d) { return y_offset_observed + y_height_observed - (Math.round(Math.log(d.observed_num) + 1)*14) ; })
    //           .on('click', function(d) {
    //              var text = "";
    //              for (k in d.observed) {
    //                 text += k + ": " + d.observed[k] + "\n";
    //              }
    //              alert(text);
    //              // alert('clicked on a cell ', d.pos);
    //           })
    //           .on('mouseover', tip.show)
    //           .on('mouseout', tip.hide)
    //           .attr("fill", "#D3FA59");
    //           // .attr("fill", "gray");
    //           // .attr("fill", "#FF9900");

    ///////////////////////////// EXPECTED DNA MUTABILITY ///////////////////////////////////////////////////////////
    var na_line = d3.svg.line()
      // .interpolate("monotone")  // basis-open
      .interpolate("step-after")  // basis-open

      .x(function(d){
          // console.log("x", x_offset + d.seq_num * (dna_width + spacer) + Math.ceil(dna_width/2));
          return x_offset + d.seq_num * (dna_width + spacer); // + Math.ceil(dna_width/2);
      })
      .y(function(d){
          // console.log("y", na_line_scale(d.mutability));
          return y_offset_dna + na_line_scale(d.mutability);
      });

    svg.append("path")
          .datum(json.na_records)
          .attr("class", "line")
          .attr("stroke", "green")
          .attr("fill", "none")
          // .attr("shape-rendering", "crispEdges")
          .attr("stroke-width", 2)
          .attr("d", na_line);


    /////////////////////////////////  DNA SEQUENCE TEXT ///////////////////
    var dnaSequence = svg.append("g")
        .attr("class", "dna-sequence")
        .selectAll(".dna-sequence")
        .data(json.na_records)
        .enter().append('text')
        .attr('y', function(d, i) {
            return y_offset_dna + y_height_dna; 
        })
        .attr('x', function(d, i) {
            return x_offset + d.seq_num * (dna_width + spacer) + dna_width / 2;
        })
        // .attr('class','matrix-label') 
        // .attr('title', function(d) {return d;})
        // .attr('font-size', '11px')
        // .style('text-anchor','middle')
        .text(function(d, i) {
            return d.seq;
        });



    
    //////////////////////////  OBSERVED DNA MUTATIONS circles //////////////////////////////////////////////////////////////
    var observed_max = d3.max(json.na_records, function(d) { return parseFloat(d.observed_num); });
    var stem_height_log_scale = d3.scale.log().base(Math.E)
        .domain([1, observed_max])
        .range([5, y_height_observed - 5]);

    svg.selectAll("stems")
              .data(json.na_records.filter(function(d){ return d.observed_num > 0;}))
          .enter().append("svg:rect")
              // .attr("class", "line")
              .attr("stroke", "#333")
              .attr("stroke-width", 1.5)
              .attr("height", function(d) {
                    return stem_height_log_scale(d.observed_num);
                 })
              .attr("width", 1)
              .attr("x", function(d) { return d3.format(".2f")(x_offset + d.seq_num * (dna_width + spacer) + Math.floor(dna_width / 2)); })
              .attr("y", function(d) { return d3.format(".2f")(y_offset_observed + y_height_observed - stem_height_log_scale(d.observed_num)); })
              .on('click', function(d) {
                 var text = "";
                 for (k in d.observed) {
                    text += k + ": " + d.observed[k] + "\n";
                 }
                 alert(text);
                 // alert('clicked on a cell ', d.pos);
              })
              .on('mouseover', tip.show)
              .on('mouseout', tip.hide);
              // .attr("fill", "#D3FA59");
              // .attr("fill", "gray");
              // .attr("fill", "#FF9900");

    svg.selectAll("circles")
              .data(json.na_records.filter(function(d){ return d.observed_num > 0;}))
          .enter().append("path")
              .attr("transform", function(d, i){
                    var x = x_offset + d.seq_num * (dna_width + spacer) + Math.ceil(dna_width / 2);
                    var y = y_offset_observed + y_height_observed - stem_height_log_scale(d.observed_num);
                    x = d3.format(".2f")(x);
                    y = d3.format(".2f")(y);
                    return "translate(" + x + "," + y + ")";
              })
              .on('click', function(d) {
                 var text = "";
                 for (k in d.observed) {
                    text += k + ": " + d.observed[k] + "\n";
                 }
                 alert(text);
                 // alert('clicked on a cell ', d.pos);
              })
              .on('mouseover', tip.show)
              .on('mouseout', tip.hide)
              .attr("d", d3.svg.symbol().type("circle").size(60))
              .attr('class', 'o-circle');



    //////////////////////////  OBSERVED DNA MUTATIONS CIRCLES //////////////////////////////////////////////////////////////
    // svg.selectAll("dot")
    //           .data(json.na_records)
    //       .enter().append("circle")
    //       .filter(function(d){ return d.observed_num > 0;})
    //           // .attr("class", "line")
    //           .attr("stroke", "gray")
    //           .attr("stroke-width", 1)
    //           .attr("r", function(d) {
    //                 if (d.observed_num == 1) {
    //                   return 5;
    //                 } else {
    //                   return 9;
    //                 }
    //              })
    //           .attr("cx", function(d) { return x_offset + d.seq_num * (dna_width + spacer) + dna_width/2; })
    //           .attr("cy", function(d) { return y_offset_dna + (y_offset_aa - y_offset_dna) / 2; })
    //           .on('click', function(d) {
    //              var text = "";
    //              for (k in d.observed) {
    //                 text += k + ": " + d.observed[k] + "\n";
    //              }
    //              alert(text);
    //              // alert('clicked on a cell ', d.pos);
    //           })
    //           .on('mouseover', tip.show)
    //           .on('mouseout', tip.hide)
    //           .attr("fill", "#B3DA29");
    //           // .attr("fill", "#FF9900");

    ////////////////////////////////////////////////////////////////////////////////////////
    svg.selectAll("triangle")
              .data(json.na_matrix.filter(function(d){ return d.observed_num > 0;}))
          .enter().append("path")
              .attr("transform", function(d, i){
                    var x = x_offset + d.x * (dna_width + spacer) + dna_width/2;
                    var y = y_offset_na_matrix + d.y * (matrix_cell_height + spacer) + (matrix_cell_height)/2;
                    return "translate(" + x + "," + y + ")";
              })
              .on('click', function(d) {
                 var text = "";
                 for (k in d.observed) {
                    text += k + ": " + d.observed[k] + "\n";
                 }
                 alert(text);
                 // alert('clicked on a cell ', d.pos);
              })
              .on('mouseover', tip.show)
              .on('mouseout', tip.hide)
              .attr("d", d3.svg.symbol().type("circle").size(20))
              .attr('class', 'o-circle');

              // .attr("stroke", "#BCBCBC")
              // .attr("fill", "#424242");
    ////////////////////////////////////////////////////////////////////////////////////////
    svg.selectAll("triangle")
              .data(json.aa_matrix.filter(function(d){ return d.observed_num > 0;}))
          .enter().append("path")
              .attr("transform", function(d, i){
                    var x = x_offset + d.x * (aa_width + spacer) + aa_width/2;
                    var y = y_offset_matrix + d.y * (matrix_cell_height + spacer) + (matrix_cell_height)/2;
                    return "translate(" + x + "," + y + ")";
              })
              .on('click', function(d) {
                 var text = "";
                 for (k in d.observed) {
                    text += k + ": " + d.observed[k] + "\n";
                 }
                 alert(text);
                 // alert('clicked on a cell ', d.pos);
              })
              .on('mouseover', tip.show)
              .on('mouseout', tip.hide)
              .attr("d", d3.svg.symbol().type("circle").size(60))
              .attr('class', 'o-circle');

              // .attr("stroke", "#BCBCBC")
              // .attr("fill", "#424242");
    ////////////////////////////////////////////////////////////////////////////////////////
    var aa_line_scale = d3.scale.linear()
        .domain([0, aa_max])
        // .range([y_offset_aa - y_offset_dna - 12, 3]);
        .range([y_height_aa + 10, 3]);

    var aa_axis = d3.svg.axis()
        .scale(aa_line_scale)
        .orient("left")
        .ticks(3);
        // .tickFormat(d3.format("e"));

    svg_labels.append("g")
        .attr("class", "mutability-axis")
        .attr("transform", "translate(" + (labels_plot_width - 2) + ", "+ (y_offset_aa) +")")
        .style("font-size", '10px')
        .style("font-weight", 'normal')
        .style("font-family", 'sans-serif')
        .call(aa_axis);

    // svg.selectAll("expected_aa")
    //       .data(json.aa_records)
    //     .enter().append("path")
    //       .attr("transform", function(d) {
    //             var x = Math.round(x_offset + d.seq_num * (aa_width + spacer));
    //             var y = Math.round(0 + aa_line_scale(d.mutability));
    //             return "translate(" + x + "," + y + ")";
    //       })
    //       // .attr("d", d3.svg.symbol().type("triangle-up"))
    //       .attr("d", "M0 0 L" + aa_width + " 0 Z")
    //       .attr("stroke-width", 2)
    //       .attr("stroke", "orange")
    //       .attr("fill", "#424242");

    // AA LINE
    var aa_line = d3.svg.line()
      .interpolate("step-after")
      .x(function(d){
          // console.log("x", x_offset + d.seq_num * (aa_width + spacer) + Math.ceil(aa_width/2));
          return x_offset + d.seq_num * (aa_width + spacer);
      })
      .y(function(d){
          // console.log("y", aa_line_scale(d.mutability));
          // return y_offset_aa + aa_line_scale(d.mutability);
          return y_offset_aa  + aa_line_scale(d.mutability);
      });

    svg.append("path")
          .datum(json.aa_records)
          .attr("class", "line")
          .attr("stroke", "orange")
          .attr("fill", "none")
          // .attr("shape-rendering", "crispEdges")
          .attr("stroke-width", 2)
          .attr("d", aa_line);


    // //label rows
    // var rowLabels = svg.append("g")
    //     .attr("class", "rowLabels")
    //     .selectAll(".rowLabel")
    //     .data(rows)
    //     .enter().append('text')
    //     .attr('y', function(d, i) {
    //         return YOffset + 2 + i*(h + spacer) + (h+spacer)/2;
    //     })
    //     .attr('x', XOffset - 5)
    //     .attr('class','matrix-label') 
    //     .attr('title', function(d) {return d;})
    //     .style('text-anchor','end')
    //     .text(function(d, i) {
    //         return d;
    //     });

    // //label columns
    // var colLabels = svg.append("g")
    //     .attr("class", "colLabels")
    //     .selectAll(".colLabel")
    //     .data(cols)
    //     .enter().append('text')
    //     .attr('x', 0)
    //     .attr('y', 0)
    //     .attr("transform", function(d, j) {
    //         return "translate("+ (XOffset + 2 + j*(w + spacer) + (w+spacer)/2) +", "+YOffset+") rotate(-45)"
    //     })
    //     .attr('class','matrix-label')
    //     .style('text-anchor','start')
    //     .text(function(d, i) {
    //         return d;
    //     });
}

function updateScanResults(signature, gene, sequence, mutations, mutation_rate, mutation_rate_value){
    $.post(
      url_api_gene, {
        'signature': signature, 'gene': gene, 'sequence': sequence, 'mutations': mutations, 'mutation_rate': mutation_rate, 'mutation_rate_value': mutation_rate_value},
        function(json) {
              $("#btn_submit").show();
              $("#btn_submit_loading").hide();
              drawScanSeqResults(json);
      }, 'json');
}


function updateScanScatterPlot(signature, gene, sequence, mutations, mutation_rate, mutation_rate_value){
    $.post(
      url_api_gene, {
        'signature': signature, 'gene': gene, 'sequence': sequence, 'mutations': mutations, 'download': '4', 'mutation_rate': mutation_rate, 'mutation_rate_value': mutation_rate_value}, // download=4 for scatter plot
        function(json) {
              $("#btn_submit").show();
              $("#btn_submit_loading").hide();
              drawScanScatterPlot(json); // scatter plot
      }, 'json');
}


function scrollToAnchor(aid){
    var aTag = $("a[name='"+ aid +"']");
    $('html,body').animate({scrollTop: aTag.offset().top},'slow');
}


function populate_list(signatures) {
    var $select = $('select[name=signature]');
    var previous_text = $select.children('option').filter(":selected").text();
    $select.find('option').remove();
    var listitems = "";
    $.each(signatures, function(key, value) {
        listitems += '<option value="' + value.id + '">' + value.name  + "</option>\n";
    });
    $select.append(listitems);
    $select.children('option').filter(function(){ return $(this).text() == previous_text; }).prop('selected', true);

    $("#download_input").val("0");

    if ($("#sequence").val() != "") {
        $("#gene_form").submit();
    }
}


function ajax_reload_list() {
    $.getJSON(uri_get_signature_list, {

              mut_type: $('input[name="sel_mut_type"]:checked').val(),
              signature_type: $('input[name="sel_signature_type"]:checked').val(),
              normalization: $('input[name="sel_normalization"]:checked').val(),

          }, function(json) {
              populate_list(json.signatures);
          });
}

var genename = "";


$(function(){

    $("#gene").bind("keydown", function (e) {
        if (e.keyCode == 13) {
            e.preventDefault();
            // console.log(e);
            $("#btn_load_seq").click();
        }
    });

    $("#seq_example1").on("click", function(){
        $("#gene").val("TP53");
        $("#btn_load_seq").click();
    });

    $("#seq_example2").on("click", function(){
        $("#gene").val("EGFR");
        $("#btn_load_seq").click();
    });

    $("#btn_load_seq").on("click", function(){
        genename = $("#gene").val().toUpperCase();
        genename = genename.replace(/[^\w\s]/gi, '')
        $("#gene").val(genename);

        $("#gene_alert").hide();
        $("#sequence").val("");
        $('#transcripts').empty();
        
        // reset the results
        resetScanSeqResults();

        if (genename === "") {
          return;
        }
        $.get(encodeURI(url_api_seq.replace('gene_var', genename)), function(data) {
            if (data.hasOwnProperty('error')) {
              $("#gene_alert_message").text(data.error);
              $("#gene_alert").toggleClass('alert-danger', true);
              $("#gene_alert").toggleClass('alert-success', false);
              $("#gene_alert").show();
            }

            if (data.hasOwnProperty('transcripts')) {
                $('#transcripts').append('<span>Transcripts: </span>');
                for (t in data.transcripts) {
                    let transcript = data.transcripts[t];
                    $('<a>',{
                        text: data.transcripts[t],
                        title: 'Load this gene transcript',
                        href: '#',
                        click: function(){ 
                          $("#gene").val(transcript);
                          $("#btn_load_seq").click();
                          return false;
                        }
                    }).appendTo('#transcripts');
                    $('#transcripts').append('<span> </span>');
                }
                $('#transcripts').append('<br/><br/>');
            }

            if (data.hasOwnProperty('sequence')) {
              // $("#sequence").attr("value", data.sequence);
              $("#sequence").val(data.sequence);
              $("#gene_alert").toggleClass('alert-danger', false);
              $("#gene_alert").toggleClass('alert-success', true);
              $("#gene_alert_message").text("Sequence for gene " + genename + " loaded");
              $("#gene_alert").show();

              $("#gene_alert").fadeTo(1000, 1).fadeOut(600);

              // analyze the results
              $("#download_input").val("0");
              $("#gene_form").submit();
            }
        });
    });


    $("#btn_submit_download_na").on("click", function(){
        $("#download_input").val("1");
        $("#gene_form").submit();
        return false;
    });

    $("#btn_submit_download_aa").on("click", function(){
        $("#download_input").val("2");
        $("#gene_form").submit();
        return false;
    });

    $("#btn_submit_download_aa_mut").on("click", function(){
        $("#download_input").val("3");
        $("#gene_form").submit();
        return false;
    });

    $("#btn_submit_show_scatterplot").on("click", function(){
        $("#download_input").val("4");
        $("#gene_form").submit();
        return false;
    });

    $("#btn_submit_download_na_mut").on("click", function(){
        $("#download_input").val("5");
        $("#gene_form").submit();
        return false;
    });

    $("#btn_submit").on("click", function(){
        $("#download_input").val("0");
        $("#gene_form").submit();
        return false;
    });

    $("#gene_form").submit(function(event){
        // console.log("DLD", $("#download_input").val());
        genename = $("#gene").val().toUpperCase();
        $("#gene").val(genename);

        if ($("#signature").val() == "") {
            alert("Signature not selected!");
            event.preventDefault();
            return;
        }

        if ($("#sequence").val() == "") {
            alert("Sequence is required!");
            event.preventDefault();
            return;
        }

        if ($("#download_input").val() == "0") {
            // force no POST submit
            event.preventDefault();
            // resetScanSeqResults();
            $("#btn_submit").hide();
            $("#btn_submit_loading").show();
            updateScanResults($("#signature").val(), $("#gene").val(), $("#sequence").val(), $("#mutations").val(), $("#mutation_rate").val(), $("#mutation_rate_value").val());
        }

        if ($("#download_input").val() == "4") {
            // force no POST submit
            event.preventDefault();
            // resetScanSeqResults();
            // $("#btn_submit").hide();
            // $("#btn_submit_loading").show();
            updateScanScatterPlot($("#signature").val(), $("#gene").val(), $("#sequence").val(), $("#mutations").val(), $("#mutation_rate").val(), $("#mutation_rate_value").val());
        }
        // // reset the status of download option
        // $("#download_input").val("0");
    });


    $("#reset_button").on('click', function(){
        resetScanSeqResults();
        $("#gene_alert").hide();

        $("#download_input").val("0");
        $("#sequence").val("");
        
        $("#mutation_rate").val("0");
        $("#mutation_rate_value").parent().hide();

        $('form')[0].reset();
    });


    $("#signature").on("change", function(){
          // analyze the results
          $("#download_input").val("0");

          if ($("#sequence").val() != "") {
              $("#gene_form").submit();
          }
    });

    $("#mutations").on("change", function(){
          // analyze the results
          $("#download_input").val("0");

          if ($("#sequence").val() != "") {
              $("#gene_form").submit();
          }
    });

    $("#mutation_rate").on("change", function(){
        $("#mutation_rate_value").parent().toggle($("#mutation_rate").val() == "2");
    });
    $("#mutation_rate_value").parent().toggle($("#mutation_rate").val() == "2");

    $("#view_mutability_button").on("click", function() {
        var accession = $("#signature").val();
        // var label = $("#signature").val();
        var label = $("#signature option:selected").text();
        // show_signature_modal(accession, label);
        show_mutability_modal(accession, label, label);
    });

    $("#view_signature_button").on("click", function() {
        var accession = $("#signature").val();
        // var label = $("#signature").val();
        var label = $("#signature option:selected").text();
        show_signature_modal(accession, label);
        // show_mutability_modal(accession, label, label);
    });

    $("input:radio").on("change", function(ev){
        var signature = $('input[name="sel_signature_type"]:checked').val();
        var selectors = $("input:radio[name=sel_mut_type]");
        if (signature == "BEN" || signature == "IS") {
            selectors.filter("[value=A]").click();
            selectors.not(":checked").parent().toggleClass("disabled", true);
            selectors.not(":checked").prop("disabled", true);
        } else {
            selectors.prop("disabled", false);
            selectors.parent().toggleClass("disabled", false);
        }
        ajax_reload_list();
    });

    // initial load
    ajax_reload_list();
    $.ajaxSetup({ cache: false });
    $("#gene_alert").hide();
    $(".results_legend").hide();
    resetScanSeqResults();
    $("#btn_submit_loading").hide();
});
