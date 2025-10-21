# Neural_Inference
NI: A light-weight header-mostly Neural Network framework for on-device Inference.

## History  
The original purpose of this framework was to provide inference for neural networks on mobile devices at at a time when there were close to none. It was implemented from scratch in C++ relying on Eigen's Library. It is completely cross-platform and currently works on Linux/Mac/Windows/iOS/Android.

## Philosophy
The idea is to keep the simplicity of the original approach while keeping efficiency and low-latency at its core.
- Fast compilation
- Mathematically easy to read
- Completely cross-platform
- Didactic: to help people learn all the 'nuts and bolts' of ML, inference, mobile devices, plugin interfaces and framework creation

## Currently working on
I want to revamp this old framework and make it truly open-source for today's rich ML ecosystem: 
- Removing dependencies to UCL copyrighted routines
- Replicating digital signal processing utilities available to current frameworks/toolkits (PyTorch, TensorFlow and Kaldi)
- Adding extensions to provide interface with specialised speech/linguistic software, such as Praat
- Adding more architectures and ops to match the current appetite for transformers
- Using flat buffers for serialisation 
- Adding web as a platform
