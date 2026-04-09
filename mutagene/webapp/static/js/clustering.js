
////////////////////////////////////////////////////////////////////////////////////////////
drawClusteringPlot = function(ref, data, points) {
    // var data = json.data;
    // console.log(data);
    // array of dict
    // string int cluster, id
    // float x, y

    var dot_size = 3;
    var width = 500;
    var height = 430;
    var margin = 80;

    d3.selectAll(ref + " > *").remove();

    var svg = d3.select(ref)
        .attr("width", width + 2 * margin)
        .attr("height", height + 15)
      .append("g")
        .attr("transform", "translate(" + (margin) + "," + 20 + ")");

    // console.log(svg);
    // d3.select("#hist2")
    //   .attr("width", 1)
    //   .attr("height", 1);

    var xValue = function(d) { return d.x;}, // data -> value
        xScale = d3.scale.linear().range([0, width]), // value -> display
        xMap = function(d) { return xScale(xValue(d));}, // data -> display
        xAxis = d3.svg.axis().scale(xScale).orient("bottom");

    var yValue = function(d) { return d.y;}, // data -> value
        yScale = d3.scale.linear().range([height, 0]), // value -> display
        yMap = function(d) { return yScale(yValue(d));}, // data -> display
        yAxis = d3.svg.axis().scale(yScale).orient("left");

    var clusters = d3.set(data.map(function(d) { return 1 + parseInt(d.cluster); })).values().sort();

    //var color = d3.scale.category10().domain([1,2,3,4,5,8,7,6,9]);
    var color = d3.scale.ordinal()
        .domain([1,2,3,4,5])
        .range(["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854"])
    // var colors = clusters.map(color);

    var min_val = Math.min(d3.min(data, yValue), d3.min(data, xValue)) - 0.01;
    var max_val = Math.max(d3.max(data, yValue), d3.max(data, xValue)) + 0.01;
    xScale.domain([min_val, max_val]);
    yScale.domain([min_val, max_val]);
    // xScale.domain([d3.min(data, xValue) - 0.0001, d3.max(data, xValue)+ 0.01] );
    // yScale.domain([d3.min(data, yValue) - 0.0001, d3.max(data, yValue)+ 0.01] );

    // x-axis
    // svg.append("g")
    //     .attr("class", "invisible-axis")
    //     .attr("transform", "translate(0," + height + ")")
    //     .call(xAxis);
      // .append("text")
      //   .attr("class", "scatterplot-label")
      //   .attr("x", width)
      //   .attr("y", 35)
      //   .style("text-anchor", "end")
      //   .text(json.xlabel);

    // y-axis
    // svg.append("g")
    //     .attr("class", "invisible-axis")
    //     .call(yAxis);
      // .append("text")
      //   .attr("class", "scatterplot-label")
      //   .attr("transform", "rotate(-90)")
      //   .attr("y", -55)
      //   .attr("dy", ".71em")
      //   // .attr("dy", "1em")
      //   .style("text-anchor", "end")
      //   .text(json.ylabel);

    // var tooltip = d3.select(".scatterplot-tooltip");
    var tip = d3.tip()
      .attr('class', 'd3-tip')
      .offset([-9, 2])
      .html(function(d) {
         var src = url_fingerprint.replace('0', d.id);
         // return 'Cluster #' + (parseInt(d.cluster) + 1) + ', sample signature ' + (d.id) + ':<br><img height="45px" width="224px" style="padding: 10px;" src="' + src + '"</img>';
         if (d.annotation) {
            return d.annotation; //+ ':<br><img height="45px" width="224px" style="padding: 10px;" src="' + src + '"</img>';
         }
         return 'Sample in cluster ' + (parseInt(d.cluster) + 1) + ':<br><img height="45px" width="224px" style="padding: 10px;" src="' + src + '"</img>';
      });
    svg.call(tip);

    svg.selectAll(".dot")
         .data(data)
           .enter().append("circle")
             .attr('class', 'cluster-dot')
             .attr('r', function(d) { return Math.min(0.8 * dot_size * Math.sqrt(d.size), 30.0);})
             .attr('cx', xMap)
             .attr('cy', yMap)
             .style('fill', function(d) {
                return color(1 + parseInt(d.cluster));
             })
             .on('mouseover', function(d){
                  // tip.show(d, document.getElementById("heading1"));
                  // console.log(this);
                  d3.select(this).classed("cluster-dot-hover", true);
                  tip.show(d, this);
             })
             .on('mouseout', function(d) {
                  d3.select(".cluster-dot-hover").classed("cluster-dot-hover", false);
                  tip.hide();
              });

    if (points) {
        svg.selectAll(".points")
            .data(points)
              .enter().append("rect")
                .attr('x', xMap)
                .attr('y', yMap)
                .attr('height', 20)
                .attr('width', 20)
                .attr('stroke', 'yellow')
                .attr('fill', 'black')
                .on('mouseover', function(d){
                     // tip.show(d, document.getElementById("heading1"));
                     // console.log(this);
                     d3.select(this).classed("cluster-dot-hover", true);
                     tip.show(d, this);
                })
                .on('mouseout', function(d) {
                     d3.select(".cluster-dot-hover").classed("cluster-dot-hover", false);
                     tip.hide();
                 });    
    }
}

////////////////////////////////////////////////////////////////////////////////////////////