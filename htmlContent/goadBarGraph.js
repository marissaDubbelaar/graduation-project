function createBarGraph(){
	// The size of the bargraph is defined.
	var margin = {top: 20, right: 20, bottom: 30, left: 40},
		width = ($("#QE_content").width() * 0.85) - margin.left - margin.right,
	    height = 500 - margin.top - margin.bottom;

	var x0 = d3.scale.ordinal()
	    .rangeRoundBands([0, width], .1);

	var x1 = d3.scale.ordinal();

	var y = d3.scale.linear()
	    .range([height, 0]);

	// Colors for the columns is defined.
	var color = d3.scale.ordinal()
	    .range(["#0080ff", "#999999", "#ffcc00"]);

	// X axis is defined.
	var xAxis = d3.svg.axis()
	    .scale(x0)
	    .orient("bottom");

	// Y axis is defined.
	var yAxis = d3.svg.axis()
	    .scale(y)
	    .orient("left")
	    .tickFormat(d3.format(".2s"));

	// The tooltip is defined.
	var tip = d3.tip()
	  .attr('class', 'd3-tip')
	  .offset([-10, 0])
	  .html(function(d) {
	    return "<strong>" + d.name + " value : </strong>" + d.value
	  })

	// The svg is drawn and added into the div with the id "TPMdiv"
	var svg = d3.select("#TPMdiv").append("svg")
	    .attr("width", width + margin.left + margin.right)
	    .attr("height", height + margin.top + margin.bottom)
	    .attr("class", "tpmValsPlot")
	  	.append("g")
	    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

	// The type of information shown (also in the legend) are defined.
	svg.call(tip);
	  var tpmTypes = ["High TPM", "Low TPM", "TPM"];

	  // The tpmVals are obtained.
	  bargraphData.forEach(function(d) {
	    d.tpmVal = tpmTypes.map(function(name) { return {name: name, value: +d[name]}; });
	  });

	  x0.domain(bargraphData.map(function(d) { return d.Gene; }));
	  x1.domain(tpmTypes).rangeRoundBands([0, x0.rangeBand()]);
	  y.domain([0, d3.max(bargraphData, function(d) { return d3.max(d.tpmVal, function(d) { return d.value; }); })]);

	  // Adds the X axis to the svg.
	  svg.append("g")
	      .attr("class", "x axis")
	      .attr("transform", "translate(0," + height + ")")
	      .call(xAxis);

	  // Adds the Y axis to the svg.
	  svg.append("g")
	      .attr("class", "y axis")
	      .call(yAxis)
	    .append("text")
	      .attr("transform", "rotate(-90)")
	      .attr("y", 6)
	      .attr("dy", ".71em")
	      .style("text-anchor", "end")
	      .text("TPM value");

	  // Defines state (used to draw all of the bars).
	  var state = svg.selectAll(".state")
	      .data(bargraphData)
	      .enter().append("g")
	      .attr("class", "state")
	      .attr("transform", function(d) { return "translate(" + x0(d.Gene) + ",0)"; });

	  // Adds each bar to the bargraph.
	  state.selectAll("rect")
	      .data(function(d) { return d.tpmVal; })
	      .enter().append("rect")
	      .attr("width", x1.rangeBand())
	      .attr("x", function(d) { return x1(d.name); })
	      .attr("y", function(d) { return y(d.value); })
	      .attr("height", function(d) { return height - y(d.value); })
	      .attr("fill", function(d) { return color(d.name); })
	      .on('mouseover', tip.show)
	      .on('mouseout', tip.hide);

	  // Adds the legend to the svg.
	  var legend = svg.selectAll(".legend")
	      .data(tpmTypes.slice().reverse())
	      .enter().append("g")
	      .attr("class", "legend")
	      .attr("transform", function(d, i) { return "translate(0," + i * 20 + ")"; });

	  // Adds the color to the bars.
	  legend.append("rect")
	      .attr("x", width - 18)
	      .attr("width", 18)
	      .attr("height", 18)
	      .style("fill", color);

	  // Adds the names on the x axis of the svg.
	  legend.append("text")
	      .attr("x", width - 24)
	      .attr("y", 9)
	      .attr("dy", ".35em")
	      .style("text-anchor", "end")
	      .text(function(d) { return d; });
}

function type(d) {
  d.tpmVal = +d.tpmVal;
  return d;
}