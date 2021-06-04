#ifndef _SZ_APS_PREPROCESSOR_HPP
#define _SZ_APS_PREPROCESSOR_HPP

#include <cstring>
#include <cassert>
#include <iostream>
#include <vector>
#include "def.hpp"
#include <bitset>
#include "preprocessor/Preprocessor.hpp"
#include "predictor/EstLorenzoPredictor.hpp"
#include "quantizer/EstQuantizer.hpp"

namespace SZ {

    template<class T, uint N>
    class APSPreprocessor : public concepts::PreprocessorInterface<T, N> {
    public:
        void process(T * data, SZ::Config<T, N>& conf) const {
        	auto quantizer = SZ::EstQuantizer<T>(conf.eb, conf.quant_bin);
        	// determine whether to transpose
        	auto global_dimensions = conf.dims;
        	size_t starting_offset = 0;
        	std::vector<size_t> offsets(conf.dims.size());
        	size_t offset = 1;
        	for(int i=global_dimensions.size() - 1; i>=0; i--){
        		offsets[i] = offset;
        		starting_offset += offset;
        		offset *= global_dimensions[i];
        		std::cout << offsets[i] << " ";
        	}
        	std::cout << std::endl;
        	std::cout << "offset = " << starting_offset << ", conf.num = " << conf.num << std::endl;

        	const int sampling_distance = 32;

        	double current_err = 0, transpose_err = 0;
        	for(int i=starting_offset; i<conf.num; i+=sampling_distance){
        		current_err += quantizer.quantize(data[i], data[i - 1]);
        		transpose_err += quantizer.quantize(data[i], data[i - offsets[0]]);
        	}
        	std::cout << "current error = " << current_err << ", transpose error = " << transpose_err << std::endl;

        	check_variance<N>(data, global_dimensions);
            // determine quantization bins

        	// determine prediction dimension
        	auto predictor = SZ::EstLorenzoPredictor3D<T>(conf.eb);
        	double est_pred_err_1d = 0;
        	double est_pred_err_2d = 0;
        	double est_pred_err_3d = 0;

        	const int block_size = 4;
            std::array<size_t, N> intra_block_dims;
            auto inter_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                    std::begin(global_dimensions),
                                                                    std::end(global_dimensions), block_size, 0);
            auto intra_block_range = std::make_shared<SZ::multi_dimensional_range<T, N>>(data,
                                                                    std::begin(global_dimensions),
                                                                    std::end(global_dimensions), 1, 0);

            {
                auto inter_begin = inter_block_range->begin();
                auto inter_end = inter_block_range->end();
                int count = 0;
                for (auto block = inter_begin; block != inter_end; ++block, ++count) {
                	// sampling every sampling_distance block
                	if(count % sampling_distance == 0){
	                    for (int i = 0; i < intra_block_dims.size(); i++) {
	                        size_t cur_index = block.get_local_index(i);
	                        size_t dims = inter_block_range->get_dimensions(i);
	                        intra_block_dims[i] = (cur_index == dims - 1) ? global_dimensions[i] - cur_index * block_size
	                                                                      : block_size;
	                    }
	                    intra_block_range->set_dimensions(intra_block_dims.begin(), intra_block_dims.end());
	                    intra_block_range->set_offsets(block.get_offset());
	                    intra_block_range->set_starting_position(block.get_local_index());

	                    auto intra_begin = intra_block_range->begin();
	                    auto intra_end = intra_block_range->end();
	                    for (auto element = intra_begin; element != intra_end; ++element) {
	                        est_pred_err_1d += fabs(*element - predictor.predict_1d(element)) + 0.5 * conf.eb;
	                        est_pred_err_2d += fabs(*element - predictor.predict_2d(element)) + 0.81 * conf.eb;
	                        est_pred_err_3d += fabs(*element - predictor.predict_3d(element)) + 1.22 * conf.eb;
	                        // est_pred_err_1d += quantizer.quantize(fabs(*element - predictor.predict_1d(element)), 0) + 0.5 * conf.eb;
	                        // est_pred_err_2d += quantizer.quantize(fabs(*element - predictor.predict_2d(element)), 0) + 0.81 * conf.eb;
	                        // est_pred_err_3d += quantizer.quantize(fabs(*element - predictor.predict_3d(element)), 0) + 1.22 * conf.eb;
	                    }
	                }
                }
            }
            std::cout << est_pred_err_1d << " " << est_pred_err_2d << " " << est_pred_err_3d << std::endl;
            if((est_pred_err_1d < est_pred_err_2d) && (est_pred_err_2d < est_pred_err_3d)){
            	conf.pred_dim = 1;
            }
            else if((est_pred_err_2d < est_pred_err_1d) && (est_pred_err_3d)){
            	conf.pred_dim = 2;
            }
            else{
            	conf.pred_dim = 3;
            }
        }
    private:
    	template<uint N1 = 3>
    	int check_variance(T * data, const std::array<size_t, 3>& dims) const {
    		const T transpose_threshold = 3.0;
    		std::vector<double> sum_x(dims[1] * dims[2], 0);
    		std::vector<double> sum_y(dims[0] * dims[2], 0);
    		std::vector<double> sum_z(dims[0] * dims[1], 0);
    		for(int i=0; i<dims[0]; i++){
    			for(int j=0; j<dims[1]; j++){
    				for(int k=0; k<dims[2]; k++){
    					T cur_data = data[i*dims[1]*dims[2] + j*dims[2] + k];
    					sum_z[i * dims[1] + j] += cur_data;
    					sum_y[i * dims[2] + k] += cur_data;
    					sum_x[j * dims[2] + k] += cur_data;
     				}
    			}
    		}
    		for(int i=0; i<dims[1] * dims[2]; i++){
    			sum_x[i] /= dims[0];
    		}
    		for(int j=0; j<dims[0] * dims[2]; j++){
    			sum_y[j] /= dims[1];
    		}
    		for(int k=0; k<dims[0] * dims[1]; k++){
    			sum_z[k] /= dims[2];
    		}
    		std::vector<double> var_x(dims[1] * dims[2], 0);
    		std::vector<double> var_y(dims[0] * dims[2], 0);
    		std::vector<double> var_z(dims[0] * dims[1], 0);
    		for(int i=0; i<dims[0]; i++){
    			for(int j=0; j<dims[1]; j++){
    				for(int k=0; k<dims[2]; k++){
    					T cur_data = data[i*dims[1]*dims[2] + j*dims[2] + k];
    					var_z[i * dims[1] + j] += (cur_data - sum_z[i * dims[1] + j])*(cur_data - sum_z[i * dims[1] + j]);
    					var_y[i * dims[2] + k] += (cur_data - sum_y[i * dims[2] + k])*(cur_data - sum_y[i * dims[2] + k]);
    					var_x[j * dims[2] + k] += (cur_data - sum_x[j * dims[2] + k])*(cur_data - sum_x[j * dims[2] + k]);
     				}
    			}
    		}
    		double std_x = 0, std_y = 0, std_z = 0;
    		for(int i=0; i<dims[1] * dims[2]; i++){
    			std_x += sqrt(var_x[i] / dims[0]);
    		}
    		std_x /= dims[1] * dims[2];
    		for(int j=0; j<dims[0] * dims[2]; j++){
    			std_y += sqrt(var_y[j] / dims[1]);
    		}
    		std_y /= dims[0] * dims[2];
    		for(int k=0; k<dims[0] * dims[1]; k++){
    			std_z += sqrt(var_z[k] / dims[2]);
    		}
    		std_z /= dims[0] * dims[1];
    		std::cout << std_x << " " << std_y << " " << std_z << std::endl;
    		// define new offsets based on standard deviation
    		size_t offsets[3];
    		if(std_y / std_x > transpose_threshold){
    			// (i, j) -> (j, i)
    			if(std_z / std_x > transpose_threshold){
    				// (i, k) -> (k, i)
    				if(std_z / std_y > transpose_threshold){
	    				// (i, j, k) -> (k, j, i)
    					offsets[0] = 1;
    					offsets[1] = dims[0];
    					offsets[2] = dims[1] * dims[0];
    					std::cout << "(k, j, i)" << std::endl;
    				}
    				else{
    					// (i, j, k) -> (j, k, i)
    					offsets[0] = 1;
    					offsets[1] = dims[2] * dims[0];
    					offsets[2] = dims[0];
    					std::cout << "(j, k, i)" << std::endl;
    				}
    			}
    			else{
    				// (i, k) -> (i, k)
					// (i, j, k) -> (j, i, k)
					offsets[0] = dims[2];
					offsets[1] = dims[2] * dims[0];
					offsets[2] = 1;
					std::cout << "(j, i, k)" << std::endl;
    			}
    		}
    		else{
    			// (i, j) -> (i, j)
				if(std_z / std_y > transpose_threshold){
    				// (j, k) -> (k, j)
    				if(std_z / std_x > transpose_threshold){
	    				// (i, j, k) -> (k, i, j)
    					offsets[0] = dims[1];
    					offsets[1] = 1;
    					offsets[2] = dims[1] * dims[0];
    					std::cout << "(k, i, j)" << std::endl;
    				}
    				else{
    					// (i, j, k) -> (i, k, j)
    					offsets[0] = dims[1] * dims[2];
    					offsets[1] = 1;
    					offsets[2] = dims[1];
    					std::cout << "(i, k, j)" << std::endl;
    				}
    			}
    			else{
    				// (j, k) -> (j, k)
    				// (i, j, k) -> (i, j, k)
    				offsets[0] = dims[1] * dims[2];
    				offsets[1] = dims[2];
    				offsets[2] = 1;
					std::cout << "(i, j, k)" << std::endl;
    			}
    		}
    		for(int i=0; i<3; i++){
    			std::cout << offsets[i] << " ";
    		}
    		std::cout << "\n";
    		// perform transpose
    		std::vector<T> tmp(dims[0] * dims[1] * dims[2]);
    		memcpy(tmp.data(), data, dims[0] * dims[1] * dims[2] * sizeof(T));
    		for(int i=0; i<dims[0]; i++){
    			for(int j=0; j<dims[1]; j++){
    				for(int k=0; k<dims[2]; k++){
    					// (i, j, k) = (std_max, std_med, std_min)
    					// data[i*dims[1]*dims[2] + j*dims[2] + k] = tmp[i*offsets[0] + j*offsets[1] + k*offsets[2]];
    					data[i*offsets[0] + j*offsets[1] + k*offsets[2]] = tmp[i*dims[1]*dims[2] + j*dims[2] + k];
     				}
    			}
    		}
    		return 0;
    	}
    };
}

#endif
