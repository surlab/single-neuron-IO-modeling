
import src.computation as comp
import numpy as np

import pytest

def my_score_function(a, b):
    return np.sum((a-b)**2)

a = np.array([1,2,3,4,5,6,7,8])

def test_fitting1_param():
	b = comp.somatic_1_fit([2], a)
	#plt.plot(a, a*1)
	#plt.scatter(a, b)
	score, m = comp.somatic_1_param(a, b, score_function=my_score_function)
	assert score < .1
	assert 1.9 < m < 2

def test_fitting2_param():
	actual_params = [1,3]
	b = comp.somatic_2_fit(actual_params, a)
	score, guessed_params = comp.somatic_2_param_lin(a, b, score_function=my_score_function)
	assert score < .1
	for i, param in enumerate(actual_params):
		assert param-.1 < guessed_params[i] <param+.1


def test_fitting_relu():
	actual_params = [1,3]
	b = comp.somatic_relu_fit(actual_params, a)
	score, guessed_params = comp.somatic_2_param_relu(a, b, score_function=my_score_function)
	assert score < .1
	for i, param in enumerate(actual_params):
		assert param-.1 < guessed_params[i] <param+.1

a2 = np.array([.1,.2,.3,.4,.5,.6,.7,.8])
def test_fitting_sigmoid():
	actual_params = [.5,.1,.02]
	b = comp.somatic_sigmoid_fit(actual_params, a2)
	score, guessed_params = comp.somatic_3_param_sigmoid(a2, b, score_function=my_score_function)
	assert score < .1
	for i, param in enumerate(actual_params):
		assert param-.1 < guessed_params[i] <param+.1

def test_normalization():
	norm_a = comp.linear_normalization(a)
	assert np.isclose(np.min(norm_a), 0)
	assert np.isclose(np.max(norm_a), 1)


def test_typical_ll():
	dummy_model_samples_uniform = [1,2]*100
	assert np.isclose(np.log(10), comp.compute_loglikelihood_discrete_distribution(dummy_model_samples_uniform, [1]))
	dummy_model_samples_uniform = [1,2]*5
	assert(
		comp.compute_loglikelihood(dummy_model_samples_uniform, [1], method='kde') == 
		comp.compute_loglikelihood(dummy_model_samples_uniform, [2], method='kde')
		)
	assert(
		comp.compute_loglikelihood(dummy_model_samples_uniform, [1.5], method='kde') >
		comp.compute_loglikelihood(dummy_model_samples_uniform, [2], method='kde')
		)
	dummy_model_samples_uniform = [1,2,3,4,5,6,7,8,9,11]*10
	assert np.isclose(np.log(1/10), comp.compute_loglikelihood_discrete_distribution(dummy_model_samples_uniform, [2]))
	assert -2.4 < comp.compute_loglikelihood(dummy_model_samples_uniform, [2], method='kde') < -2.2



def test_raises_ll():
    with pytest.raises(AssertionError):
        comp.compute_loglikelihood_discrete_distribution([1,2], [1,2,3])

    with pytest.raises(AssertionError):
        comp.compute_loglikelihood_discrete_distribution([1,2,3,4], [1,2,3])