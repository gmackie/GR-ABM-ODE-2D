/*
 * isolinesvisualization.h
 *
 *  Created on: 7-nov-2008
 *      Author: s030858
 */

#ifndef ISOLINESVISUALIZATION_H_
#define ISOLINESVISUALIZATION_H_

#include "visualization.h"
#include "vectordatasets/vector.h"
#include "scalardatasets/scalarnormalizer.h"
#include "scalardatasets/scalargrid.h"

class IsolinesVisualization: public Visualization
{
public:
	IsolinesVisualization(int DIM, const ScalarNormalizer* pScalarNormalizer, const ScalarGrid* pScalarGrid);
	virtual ~IsolinesVisualization();
	void visualize(bool blend, const Simulation* pSimulation, const ColorMap* pColorMap) const;
	const std::vector<float>& getTargetValueVector() const;
	void setTargetValueVector(const std::vector<float>& targetValueVector);
	float getLineWidth() const;
	void setLineWidth(float lineWidth);

private:
	void processTriangle(const vec2f& a, const vec2f& b, const vec2f& c,
			float val_a, float val_b, float val_c, float targetValue) const;
	vec2f interpolate(const vec2f& p_i, const vec2f& p_j, float v_i, float v_j, float v) const;
	const ScalarNormalizer* _pScalarNormalizer;
	const ScalarGrid* _pScalarGrid;
	std::vector<float> _targetValueVector;
	float _lineWidth;
};

inline float IsolinesVisualization::getLineWidth() const
{
	return _lineWidth;
}

inline void IsolinesVisualization::setLineWidth(float lineWidth)
{
	_lineWidth = lineWidth;
}

inline const std::vector<float>& IsolinesVisualization::getTargetValueVector() const
{
	return _targetValueVector;
}

inline void IsolinesVisualization::setTargetValueVector(const std::vector<float>& targetValueVector)
{
	_targetValueVector = targetValueVector;
}

#endif /* ISOLINESVISUALIZATION_H_ */
