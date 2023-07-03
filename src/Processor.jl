abstract type AbstractProcessor end

mutable struct Processor <: AbstractProcessor
    process_function::Function
end
(p::Processor)(record) = p.process_function(record)
(p::Processor)(record1, record2) = p.process_function(record1, record2)

mutable struct ExternalTool <: AbstractProcessor
    process_function::Function
end
(p::ExternalTool)(x) = p.process_function(x)
(p::ExternalTool)(x1,x2) = p.process_function(x1,x2)